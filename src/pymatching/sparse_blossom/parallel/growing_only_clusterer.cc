// Copyright 2026
// Experimental sparse-blossom-like growing-only clustering path.

#include "pymatching/sparse_blossom/parallel/growing_only_clusterer.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "pymatching/sparse_blossom/flooder/detector_node.h"

namespace pm {
namespace {
using Clock = std::chrono::steady_clock;

uint64_t elapsed_ns(Clock::time_point start) {
    const char* disable = std::getenv("PYMATCHING_DISABLE_WALL_TIMING");
    const char* event_only_disable = std::getenv("PYMATCHING_EVENTCOUNT_ONLY_NO_WALL_TIMING");
    if ((disable != nullptr && std::string(disable) == "1") ||
        (event_only_disable != nullptr && std::string(event_only_disable) == "1")) {
        return 0;
    }
    return static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count());
}

bool env_is_set_to(const char* name, const char* expected) {
    const char* value = std::getenv(name);
    return value != nullptr && std::string(value) == expected;
}

bool env_is_enabled(const char* name) {
    return env_is_set_to(name, "1");
}

bool env_is_disabled(const char* name) {
    return env_is_set_to(name, "0");
}

std::vector<uint64_t> normalise_active_detectors_for_growing_only(
    const MatchingGraph& graph,
    std::vector<uint64_t> active_detectors) {
    // The hierarchy builder keeps the residual syndrome sorted and unique.
    // Preserve that common case without allocating/copying another vector or
    // sorting again.  If a caller supplies unsorted input, fall back to the
    // old normalize-by-sort behavior.
    bool sorted_unique = true;
    bool have_last = false;
    uint64_t last = 0;
    size_t out = 0;
    for (size_t i = 0; i < active_detectors.size(); i++) {
        uint64_t detector = active_detectors[i];
        if (detector >= graph.nodes.size()) {
            throw std::invalid_argument(
                "The detection event with index " + std::to_string(detector) +
                " does not correspond to a node in the graph, which only has " +
                std::to_string(graph.nodes.size()) + " nodes.");
        }
        if (detector < graph.is_user_graph_boundary_node.size() && graph.is_user_graph_boundary_node[detector]) {
            sorted_unique = false;
            continue;
        }
        if (have_last && detector <= last) {
            sorted_unique = false;
        }
        have_last = true;
        last = detector;
        active_detectors[out++] = detector;
    }
    active_detectors.resize(out);
    if (!sorted_unique) {
        std::sort(active_detectors.begin(), active_detectors.end());
        active_detectors.erase(std::unique(active_detectors.begin(), active_detectors.end()), active_detectors.end());
    }
    return active_detectors;
}

cumulative_time_int cached_distance_from_processing_graph_cache(
    const ProcessingClusterGraphCache& graph_cache,
    size_t source,
    size_t destination) {
    if (!graph_cache.has_all_pairs_distances) {
        return GROWING_ONLY_INF_DISTANCE;
    }
    if (source > destination) {
        std::swap(source, destination);
    }
    if (source + 1 >= graph_cache.packed_all_pairs_row_offsets.size()) {
        return GROWING_ONLY_INF_DISTANCE;
    }
    uint64_t index = graph_cache.packed_all_pairs_row_offsets[source] + (destination - source);
    if (index >= graph_cache.packed_all_pairs_interior_distances.size()) {
        return GROWING_ONLY_INF_DISTANCE;
    }
    auto packed = graph_cache.packed_all_pairs_interior_distances[index];
    return packed == std::numeric_limits<uint32_t>::max() ? GROWING_ONLY_INF_DISTANCE : static_cast<cumulative_time_int>(packed);
}

cumulative_time_int min_interior_edge_weight(const MatchingGraph& graph) {
    cumulative_time_int best = GROWING_ONLY_INF_DISTANCE;
    for (const auto& node : graph.nodes) {
        for (size_t k = 0; k < node.neighbors.size(); k++) {
            if (node.neighbors[k] == nullptr) {
                continue;
            }
            best = std::min(best, static_cast<cumulative_time_int>(node.neighbor_weights[k]));
        }
    }
    return best;
}

cumulative_time_int cached_min_interior_edge_weight(const MatchingGraph& graph) {
    static thread_local std::unordered_map<const MatchingGraph*, cumulative_time_int> cache;
    auto it = cache.find(&graph);
    if (it != cache.end()) {
        return it->second;
    }
    auto value = min_interior_edge_weight(graph);
    cache.emplace(&graph, value);
    return value;
}

cumulative_time_int cached_min_nearest_boundary_match_distance(
    const ProcessingClusterGraphCache* graph_cache) {
    if (graph_cache == nullptr || graph_cache->nearest_boundary_match_distance_by_vertex.empty()) {
        return GROWING_ONLY_INF_DISTANCE;
    }
    static thread_local std::unordered_map<const ProcessingClusterGraphCache*, cumulative_time_int> cache;
    auto it = cache.find(graph_cache);
    if (it != cache.end()) {
        return it->second;
    }
    cumulative_time_int best = GROWING_ONLY_INF_DISTANCE;
    for (auto distance : graph_cache->nearest_boundary_match_distance_by_vertex) {
        best = std::min(best, distance);
    }
    cache.emplace(graph_cache, best);
    return best;
}

}  // namespace

bool GrowingOnlyClusterer::PairCollisionEvent::operator<(const PairCollisionEvent& other) const {
    if (time != other.time) return time < other.time;
    if (distance != other.distance) return distance < other.distance;
    if (a != other.a) return a < other.a;
    return b < other.b;
}

GrowingOnlyClusterer::GrowingOnlyClusterer(
    const MatchingGraph& graph,
    std::vector<uint64_t> active_detectors,
    size_t level,
    const ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache)
    : graph(graph),
      active_detectors(normalise_active_detectors_for_growing_only(graph, std::move(active_detectors))),
      level(level),
      config(config),
      graph_cache(graph_cache),
      bounds(processing_cluster_bounds_for_level(level, config)) {
    if (level == 0) {
        throw std::invalid_argument("GrowingOnlyClusterer level must be at least 1.");
    }
}

std::vector<cumulative_time_int> GrowingOnlyClusterer::dijkstra_from_source(uint64_t source) const {
    std::vector<cumulative_time_int> distances(graph.nodes.size(), GROWING_ONLY_INF_DISTANCE);
    using QueueEntry = std::pair<cumulative_time_int, size_t>;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
    distances[source] = 0;
    queue.push({0, static_cast<size_t>(source)});
    while (!queue.empty()) {
        auto [distance, node] = queue.top();
        queue.pop();
        if (distance != distances[node]) {
            continue;
        }
        const auto& detector_node = graph.nodes[node];
        for (size_t k = 0; k < detector_node.neighbors.size(); k++) {
            const DetectorNode* neighbor = detector_node.neighbors[k];
            if (neighbor == nullptr) {
                continue;
            }
            auto neighbor_index = static_cast<size_t>(neighbor - graph.nodes.data());
            auto edge_weight = static_cast<cumulative_time_int>(detector_node.neighbor_weights[k]);
            if (distance > GROWING_ONLY_INF_DISTANCE - edge_weight) {
                continue;
            }
            auto candidate = distance + edge_weight;
            if (candidate < distances[neighbor_index]) {
                distances[neighbor_index] = candidate;
                queue.push({candidate, neighbor_index});
            }
        }
    }
    return distances;
}

cumulative_time_int GrowingOnlyClusterer::distance_between_active_detectors(size_t a, size_t b) const {
    if (graph_cache != nullptr && graph_cache->has_all_pairs_distances) {
        return cached_distance_from_processing_graph_cache(
            *graph_cache,
            static_cast<size_t>(active_detectors[a]),
            static_cast<size_t>(active_detectors[b]));
    }
    return GROWING_ONLY_INF_DISTANCE;
}

std::vector<GrowingOnlyClusterer::PairCollisionEvent> GrowingOnlyClusterer::generate_pair_collision_events() {
    std::vector<PairCollisionEvent> events;
    if (active_detectors.size() < 2) {
        return events;
    }
    events.reserve((active_detectors.size() * (active_detectors.size() - 1)) / 2);

    if (graph_cache != nullptr && graph_cache->has_all_pairs_distances) {
        for (size_t i = 0; i < active_detectors.size(); i++) {
            for (size_t j = i + 1; j < active_detectors.size(); j++) {
                auto distance = distance_between_active_detectors(i, j);
                if (distance >= GROWING_ONLY_INF_DISTANCE) {
                    continue;
                }
                events.push_back(PairCollisionEvent{
                    static_cast<cumulative_time_int>(distance / 2),
                    distance,
                    static_cast<uint32_t>(i),
                    static_cast<uint32_t>(j),
                });
            }
        }
    } else {
        // Slow fallback for standalone experiments without the all-pairs graph cache.
        for (size_t i = 0; i < active_detectors.size(); i++) {
            auto distances = dijkstra_from_source(active_detectors[i]);
            for (size_t j = i + 1; j < active_detectors.size(); j++) {
                auto distance = distances[active_detectors[j]];
                if (distance >= GROWING_ONLY_INF_DISTANCE) {
                    continue;
                }
                events.push_back(PairCollisionEvent{
                    static_cast<cumulative_time_int>(distance / 2),
                    distance,
                    static_cast<uint32_t>(i),
                    static_cast<uint32_t>(j),
                });
            }
        }
    }

    std::sort(events.begin(), events.end());
    return events;
}


std::vector<GrowingOnlyClusterer::PairCollisionEvent> GrowingOnlyClusterer::generate_sparse_frontier_collision_events() {
    std::vector<PairCollisionEvent> events;
    if (active_detectors.size() < 2) {
        return events;
    }

    // Sparse-blossom-like growing-only frontier pass.  This mirrors the existing
    // lightweight threshold component flooder in processing_cluster.cc: all
    // active detectors grow simultaneously up to buffer_bound.  Whenever two
    // different source labels meet across a detector-graph edge (or one source
    // reaches another active detector), we have a collision witness.  These
    // witnesses are enough to build the same buffer-threshold components as the
    // existing preprocessing path, while avoiding a sorted all-active-pairs event
    // list.
    std::vector<cumulative_time_int> distances(graph.nodes.size(), GROWING_ONLY_INF_DISTANCE);
    std::vector<size_t> owner(graph.nodes.size(), SIZE_MAX);
    std::vector<int64_t> source_detector_index(graph.nodes.size(), -1);
    using QueueEntry = std::tuple<cumulative_time_int, size_t, size_t>;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;

    for (size_t source_index = 0; source_index < active_detectors.size(); source_index++) {
        auto source = static_cast<size_t>(active_detectors[source_index]);
        distances[source] = 0;
        owner[source] = source_index;
        source_detector_index[source] = static_cast<int64_t>(source_index);
        queue.push({0, source, source_index});
    }

    while (!queue.empty()) {
        auto [distance, node, source_index] = queue.top();
        queue.pop();
        if (distance != distances[node] || source_index != owner[node]) {
            continue;
        }
        if (distance > bounds.buffer_bound) {
            continue;
        }
        const DetectorNode& detector_node = graph.nodes[node];
        for (size_t k = 0; k < detector_node.neighbors.size(); k++) {
            const DetectorNode* neighbor = detector_node.neighbors[k];
            if (neighbor == nullptr) {
                continue;
            }
            auto neighbor_index = static_cast<size_t>(neighbor - graph.nodes.data());
            auto edge_weight = static_cast<cumulative_time_int>(detector_node.neighbor_weights[k]);
            if (distance > GROWING_ONLY_INF_DISTANCE - edge_weight) {
                continue;
            }
            auto candidate = distance + edge_weight;
            if (candidate > bounds.buffer_bound) {
                continue;
            }
            if (candidate < distances[neighbor_index] ||
                (candidate == distances[neighbor_index] && source_index < owner[neighbor_index])) {
                distances[neighbor_index] = candidate;
                owner[neighbor_index] = source_index;
                queue.push({candidate, neighbor_index, source_index});
            }
        }
    }

    events.reserve(active_detectors.size());
    auto add_event = [&](size_t a, size_t b, cumulative_time_int distance) {
        if (a == b || distance > bounds.buffer_bound || distance >= GROWING_ONLY_INF_DISTANCE) {
            return;
        }
        if (a > b) {
            std::swap(a, b);
        }
        events.push_back(PairCollisionEvent{
            static_cast<cumulative_time_int>(distance / 2),
            distance,
            static_cast<uint32_t>(a),
            static_cast<uint32_t>(b),
        });
    };

    for (size_t node = 0; node < graph.nodes.size(); node++) {
        auto owning_source = owner[node];
        if (owning_source == SIZE_MAX || distances[node] > bounds.buffer_bound) {
            continue;
        }
        auto self_source = source_detector_index[node];
        if (self_source >= 0 && static_cast<size_t>(self_source) != owning_source) {
            add_event(owning_source, static_cast<size_t>(self_source), distances[node]);
        }
        const DetectorNode& detector_node = graph.nodes[node];
        for (size_t k = 0; k < detector_node.neighbors.size(); k++) {
            const DetectorNode* neighbor = detector_node.neighbors[k];
            if (neighbor == nullptr) {
                continue;
            }
            auto neighbor_index = static_cast<size_t>(neighbor - graph.nodes.data());
            if (neighbor_index <= node) {
                continue;
            }
            auto neighbor_owner = owner[neighbor_index];
            if (neighbor_owner == SIZE_MAX || neighbor_owner == owning_source || distances[neighbor_index] > bounds.buffer_bound) {
                continue;
            }
            auto edge_weight = static_cast<cumulative_time_int>(detector_node.neighbor_weights[k]);
            if (distances[node] > GROWING_ONLY_INF_DISTANCE - edge_weight ||
                distances[node] + edge_weight > GROWING_ONLY_INF_DISTANCE - distances[neighbor_index]) {
                continue;
            }
            auto witness_distance = distances[node] + edge_weight + distances[neighbor_index];
            add_event(owning_source, neighbor_owner, witness_distance);
        }
    }

    std::sort(events.begin(), events.end());
    events.erase(
        std::unique(
            events.begin(),
            events.end(),
            [](const PairCollisionEvent& a, const PairCollisionEvent& b) {
                return a.a == b.a && a.b == b.b && a.distance == b.distance;
            }),
        events.end());
    return events;
}

bool GrowingOnlyClusterer::use_sparse_frontier_path() const {
    if (env_is_enabled("PYMATCHING_GROWING_ONLY_FULL_PAIR_EVENTS") ||
        env_is_enabled("PYMATCHING_GROWING_ONLY_REFERENCE_ACTIVE_PAIR_TABLE")) {
        return false;
    }
    if (env_is_disabled("PYMATCHING_GROWING_ONLY_SPARSE_FRONTIER")) {
        return false;
    }
    return true;
}

bool GrowingOnlyClusterer::use_reference_active_pair_table_path() const {
    return env_is_enabled("PYMATCHING_GROWING_ONLY_REFERENCE_ACTIVE_PAIR_TABLE") ||
           env_is_enabled("PYMATCHING_GROWING_ONLY_FULL_PAIR_EVENTS");
}

bool GrowingOnlyClusterer::skip_forced_max_exact_metadata() const {
    if (!force_at_max_level()) {
        return false;
    }
    if (env_is_enabled("PYMATCHING_GROWING_ONLY_FORCE_EXACT_AT_MAX_LEVEL")) {
        return false;
    }
    if (env_is_disabled("PYMATCHING_SKIP_FORCED_MAX_EXACT_CLUSTER_METADATA")) {
        return false;
    }
    return true;
}

bool GrowingOnlyClusterer::compute_exact_external_metadata() const {
    return env_is_enabled("PYMATCHING_GROWING_ONLY_EXACT_EXTERNAL_METADATA");
}

void GrowingOnlyClusterer::compute_exact_component_pair_metadata() {
    for (size_t i = 0; i < meta.size(); i++) {
        auto root = find(static_cast<uint32_t>(i));
        if (root != i || members_by_root[root].empty()) {
            continue;
        }
        meta[root].max_internal_pair_collision_time = 0;
        meta[root].max_internal_pair_distance = 0;
        internal_collision_count_by_root[root] = 0;
    }
    stats.internal_collisions_recorded = 0;
    stats.exact_internal_pair_checks = 0;
    stats.exact_diameter_checks = 0;

    const bool forced = force_at_max_level();
    if (skip_forced_max_exact_metadata()) {
        for (size_t root_index = 0; root_index < meta.size(); root_index++) {
            auto root = find(static_cast<uint32_t>(root_index));
            if (root != root_index || members_by_root[root].empty()) {
                continue;
            }
            stats.forced_max_level_shortcut_count++;
            stats.max_component_active_size = std::max(stats.max_component_active_size, members_by_root[root].size());
        }
        return;
    }

    const size_t max_active_detectors = max_active_detectors_for_this_level();
    for (size_t root_index = 0; root_index < meta.size(); root_index++) {
        auto root = find(static_cast<uint32_t>(root_index));
        if (root != root_index || members_by_root[root].empty()) {
            continue;
        }
        auto& members = members_by_root[root];
        stats.max_component_active_size = std::max(stats.max_component_active_size, members.size());

        // Candidate-only exact diameter: if a component is going to be rejected
        // for a cheap reason, do not spend O(m^2) distance checks on it.
        if (!forced && max_active_detectors != 0 && members.size() > max_active_detectors) {
            continue;
        }
        if (!forced && !parity_can_stop_locally(root)) {
            continue;
        }
        if (members.size() < 2) {
            continue;
        }

        stats.exact_diameter_checks++;
        bool exceeded_bound = false;
        for (size_t a = 0; a < members.size() && !exceeded_bound; a++) {
            std::vector<cumulative_time_int> distances;
            if (graph_cache == nullptr || !graph_cache->has_all_pairs_distances) {
                distances = dijkstra_from_source(active_detectors[members[a]]);
            }
            for (size_t b = a + 1; b < members.size(); b++) {
                cumulative_time_int distance;
                if (graph_cache != nullptr && graph_cache->has_all_pairs_distances) {
                    distance = distance_between_active_detectors(members[a], members[b]);
                } else {
                    distance = distances[active_detectors[members[b]]];
                }
                stats.exact_internal_pair_checks++;
                if (distance >= GROWING_ONLY_INF_DISTANCE) {
                    meta[root].max_internal_pair_distance = GROWING_ONLY_INF_DISTANCE;
                    meta[root].max_internal_pair_collision_time = GROWING_ONLY_INF_DISTANCE / 2;
                    exceeded_bound = !forced;
                    break;
                }
                auto t = static_cast<cumulative_time_int>(distance / 2);
                meta[root].max_internal_pair_distance = std::max(meta[root].max_internal_pair_distance, distance);
                meta[root].max_internal_pair_collision_time = std::max(meta[root].max_internal_pair_collision_time, t);
                internal_collision_count_by_root[root]++;
                stats.internal_collisions_recorded++;
                if (!forced && distance > bounds.diameter_bound) {
                    // Acceptance only needs to know that the bound is exceeded.
                    // Stop early instead of completing the dense component table.
                    exceeded_bound = true;
                    break;
                }
            }
        }
    }
}

void GrowingOnlyClusterer::compute_first_external_collision_times() {
    for (size_t i = 0; i < meta.size(); i++) {
        auto root = find(static_cast<uint32_t>(i));
        if (root != i || members_by_root[root].empty()) {
            continue;
        }
        meta[root].first_external_collision_time = GROWING_ONLY_INF_DISTANCE;
        external_collision_count_by_root[root] = 0;
    }
    stats.external_collisions_recorded = 0;
    stats.exact_external_pair_checks = 0;

    for (size_t i = 0; i < active_detectors.size(); i++) {
        auto ri = find(static_cast<uint32_t>(i));
        for (size_t j = i + 1; j < active_detectors.size(); j++) {
            auto rj = find(static_cast<uint32_t>(j));
            if (ri == rj) {
                continue;
            }
            auto distance = distance_between_active_detectors(i, j);
            if (distance >= GROWING_ONLY_INF_DISTANCE) {
                continue;
            }
            auto t = static_cast<cumulative_time_int>(distance / 2);
            if (t >= meta[ri].ready_time) {
                meta[ri].first_external_collision_time = std::min(meta[ri].first_external_collision_time, t);
            }
            if (t >= meta[rj].ready_time) {
                meta[rj].first_external_collision_time = std::min(meta[rj].first_external_collision_time, t);
            }
            external_collision_count_by_root[ri]++;
            external_collision_count_by_root[rj]++;
            stats.external_collisions_recorded++;
            stats.exact_external_pair_checks++;
        }
    }
}



bool GrowingOnlyClusterer::run_active_pair_table_component_path(GrowingOnlyClustererResult& result) {
    if (graph_cache == nullptr || !graph_cache->has_all_pairs_distances || active_detectors.empty()) {
        return false;
    }
    if (!use_reference_active_pair_table_path()) {
        return false;
    }

    auto generation_start = Clock::now();
    std::vector<PairCollisionEvent> pairs;
    pairs.reserve((active_detectors.size() * (active_detectors.size() - 1)) / 2);
    for (size_t i = 0; i < active_detectors.size(); i++) {
        for (size_t j = i + 1; j < active_detectors.size(); j++) {
            auto distance = distance_between_active_detectors(i, j);
            if (distance >= GROWING_ONLY_INF_DISTANCE) {
                continue;
            }
            pairs.push_back(PairCollisionEvent{
                static_cast<cumulative_time_int>(distance / 2),
                distance,
                static_cast<uint32_t>(i),
                static_cast<uint32_t>(j),
            });
        }
    }
    stats.event_generation_wall_ns += elapsed_ns(generation_start);
    stats.used_active_pair_table = true;
    stats.generated_pair_collision_events = 0;

    auto processing_start = Clock::now();
    std::vector<uint32_t> parent(active_detectors.size());
    std::vector<uint8_t> rank(active_detectors.size(), 0);
    std::iota(parent.begin(), parent.end(), 0);
    auto local_find = [&](uint32_t k) {
        while (parent[k] != k) {
            parent[k] = parent[parent[k]];
            k = parent[k];
        }
        return k;
    };
    auto local_unite = [&](uint32_t a, uint32_t b) {
        a = local_find(a);
        b = local_find(b);
        if (a == b) {
            return false;
        }
        if (rank[a] < rank[b]) {
            std::swap(a, b);
        }
        parent[b] = a;
        if (rank[a] == rank[b]) {
            rank[a]++;
        }
        return true;
    };
    for (const auto& pair : pairs) {
        if (pair.distance <= bounds.buffer_bound) {
            local_unite(pair.a, pair.b);
        }
    }

    std::vector<std::vector<uint32_t>> components_by_root(active_detectors.size());
    for (size_t i = 0; i < active_detectors.size(); i++) {
        components_by_root[local_find(static_cast<uint32_t>(i))].push_back(static_cast<uint32_t>(i));
    }
    std::vector<uint32_t> class_root_by_active(active_detectors.size());
    std::vector<uint32_t> local_index_by_active(active_detectors.size());
    for (size_t r = 0; r < components_by_root.size(); r++) {
        auto& component = components_by_root[r];
        if (component.empty()) {
            continue;
        }
        std::sort(component.begin(), component.end());
        uint32_t root = component[0];
        members_by_root[root] = component;
        meta[root].uf_parent = root;
        meta[root].size = static_cast<uint32_t>(component.size());
        meta[root].ready_time = 0;
        meta[root].max_internal_pair_collision_time = 0;
        meta[root].max_internal_pair_distance = 0;
        meta[root].first_external_collision_time = GROWING_ONLY_INF_DISTANCE;
        meta[root].nearest_boundary_match_distance = GROWING_ONLY_INF_DISTANCE;
        internal_collision_count_by_root[root] = 0;
        external_collision_count_by_root[root] = 0;
        for (size_t k = 0; k < component.size(); k++) {
            uint32_t member = component[k];
            meta[member].uf_parent = root;
            class_root_by_active[member] = root;
            local_index_by_active[member] = static_cast<uint32_t>(k);
            if (member != root) {
                members_by_root[member].clear();
                stats.union_count++;
            }
            meta[root].nearest_boundary_match_distance = std::min(
                meta[root].nearest_boundary_match_distance,
                meta[member].nearest_boundary_match_distance);
        }
    }

    std::vector<std::vector<PairCollisionEvent>> internal_edges_by_root(active_detectors.size());
    stats.internal_collisions_recorded = 0;
    stats.external_collisions_recorded = 0;
    stats.exact_internal_pair_checks = 0;
    stats.exact_external_pair_checks = 0;
    for (const auto& pair : pairs) {
        uint32_t ra = class_root_by_active[pair.a];
        uint32_t rb = class_root_by_active[pair.b];
        if (ra == rb) {
            meta[ra].max_internal_pair_distance = std::max(meta[ra].max_internal_pair_distance, pair.distance);
            meta[ra].max_internal_pair_collision_time = std::max(meta[ra].max_internal_pair_collision_time, pair.time);
            internal_collision_count_by_root[ra]++;
            stats.internal_collisions_recorded++;
            stats.exact_internal_pair_checks++;
            if (pair.distance <= bounds.buffer_bound) {
                internal_edges_by_root[ra].push_back(PairCollisionEvent{
                    pair.time,
                    pair.distance,
                    local_index_by_active[pair.a],
                    local_index_by_active[pair.b],
                });
            }
        } else {
            if (pair.time >= meta[ra].ready_time) {
                meta[ra].first_external_collision_time = std::min(meta[ra].first_external_collision_time, pair.time);
            }
            if (pair.time >= meta[rb].ready_time) {
                meta[rb].first_external_collision_time = std::min(meta[rb].first_external_collision_time, pair.time);
            }
            external_collision_count_by_root[ra]++;
            external_collision_count_by_root[rb]++;
            stats.external_collisions_recorded++;
            stats.exact_external_pair_checks++;
        }
    }

    for (size_t r = 0; r < internal_edges_by_root.size(); r++) {
        if (members_by_root[r].empty() || internal_edges_by_root[r].empty()) {
            continue;
        }
        auto& edges = internal_edges_by_root[r];
        std::sort(edges.begin(), edges.end());
        std::vector<uint32_t> mst_parent(members_by_root[r].size());
        std::iota(mst_parent.begin(), mst_parent.end(), 0);
        auto mst_find = [&](uint32_t k) {
            while (mst_parent[k] != k) {
                mst_parent[k] = mst_parent[mst_parent[k]];
                k = mst_parent[k];
            }
            return k;
        };
        size_t connected_edges = 0;
        for (const auto& edge : edges) {
            uint32_t a = mst_find(edge.a);
            uint32_t b = mst_find(edge.b);
            if (a == b) {
                continue;
            }
            mst_parent[b] = a;
            meta[r].ready_time = std::max(meta[r].ready_time, edge.time);
            if (++connected_edges + 1 == members_by_root[r].size()) {
                break;
            }
        }
    }

    stats.processed_pair_collision_events = pairs.size();
    evaluate_final_components(result);
    stats.event_processing_wall_ns += elapsed_ns(processing_start);
    result.stats = stats;
    collect_residual(result);
    return true;
}


bool GrowingOnlyClusterer::run_forced_max_level_direct_component_lookup_path(GrowingOnlyClustererResult& result) {
    if (!force_at_max_level() || !skip_forced_max_exact_metadata() || graph_cache == nullptr || active_detectors.empty()) {
        return false;
    }

    const bool buffer_covers_all_finite_pairs =
        graph_cache->has_all_pairs_distances && graph_cache->all_interior_pairs_finite &&
        graph_cache->max_finite_interior_distance < GROWING_ONLY_INF_DISTANCE &&
        bounds.buffer_bound >= graph_cache->max_finite_interior_distance;
    const bool can_group_by_interior_component =
        graph_cache->interior_component_id_by_vertex.size() == graph_cache->num_nodes;
    if (!buffer_covers_all_finite_pairs && !can_group_by_interior_component) {
        return false;
    }

    auto start = Clock::now();
    stats.used_component_lookup = true;
    stats.max_residual_active_size = std::max(stats.max_residual_active_size, active_detectors.size());

    auto emit_direct_cluster = [&](std::vector<uint64_t>&& detectors) {
        if (detectors.empty()) {
            return;
        }
        GrowingOnlyAcceptedCluster cluster;
        cluster.active_detectors = std::move(detectors);
        // active_detectors is already normalized/sorted; grouped slices are sorted
        // because they are produced in detector order.
        cluster.ready_time = 0;
        cluster.max_internal_pair_collision_time = 0;
        cluster.max_internal_pair_distance = 0;
        cluster.first_external_collision_time = GROWING_ONLY_INF_DISTANCE;
        cluster.nearest_boundary_match_distance = GROWING_ONLY_INF_DISTANCE;
        cluster.parity_can_stop_locally = true;
        cluster.forced_by_max_level = true;
        cluster.diameter = 0;
        cluster.buffer = GROWING_ONLY_INF_DISTANCE;
        stats.max_component_active_size = std::max(stats.max_component_active_size, cluster.active_detectors.size());
        stats.forced_max_level_shortcut_count++;
        stats.forced_by_max_level_clusters++;
        stats.removed_or_hidden_clusters++;
        result.accepted_clusters.push_back(std::move(cluster));
    };

    if (buffer_covers_all_finite_pairs) {
        emit_direct_cluster(std::move(active_detectors));
        stats.event_generation_wall_ns += elapsed_ns(start);
        result.stats = stats;
        return true;
    }

    // Fast common case for connected detector graphs: all residual detections
    // share the same interior component id. Avoid constructing per-detector UF
    // metadata and avoid an unordered_map in the forced max-level path.
    size_t first_component_id = SIZE_MAX;
    bool all_same_component = true;
    for (size_t i = 0; i < active_detectors.size(); i++) {
        size_t detector = static_cast<size_t>(active_detectors[i]);
        size_t component_id = detector < graph_cache->interior_component_id_by_vertex.size()
                                  ? graph_cache->interior_component_id_by_vertex[detector]
                                  : SIZE_MAX;
        if (i == 0) {
            first_component_id = component_id;
        } else if (component_id != first_component_id) {
            all_same_component = false;
            break;
        }
    }
    if (all_same_component && first_component_id != SIZE_MAX) {
        emit_direct_cluster(std::move(active_detectors));
        stats.event_generation_wall_ns += elapsed_ns(start);
        result.stats = stats;
        return true;
    }

    // The common rotated-code graph has only a small number of interior
    // connected components.  Group in one pass instead of materializing and
    // sorting (component_id, detector) pairs on every max-level shot.
    std::vector<size_t> component_ids;
    std::vector<std::vector<uint64_t>> grouped_detectors;
    component_ids.reserve(4);
    grouped_detectors.reserve(4);
    for (auto detector_u64 : active_detectors) {
        size_t detector = static_cast<size_t>(detector_u64);
        size_t component_id = detector < graph_cache->interior_component_id_by_vertex.size()
                                  ? graph_cache->interior_component_id_by_vertex[detector]
                                  : SIZE_MAX;
        if (component_id == SIZE_MAX) {
            std::vector<uint64_t> singleton;
            singleton.push_back(detector_u64);
            emit_direct_cluster(std::move(singleton));
            continue;
        }
        size_t out = SIZE_MAX;
        for (size_t k = 0; k < component_ids.size(); k++) {
            if (component_ids[k] == component_id) {
                out = k;
                break;
            }
        }
        if (out == SIZE_MAX) {
            component_ids.push_back(component_id);
            grouped_detectors.emplace_back();
            out = grouped_detectors.size() - 1;
        }
        grouped_detectors[out].push_back(detector_u64);
    }
    for (auto& component : grouped_detectors) {
        emit_direct_cluster(std::move(component));
    }

    active_detectors.clear();
    stats.event_generation_wall_ns += elapsed_ns(start);
    result.stats = stats;
    return true;
}

bool GrowingOnlyClusterer::run_precomputed_component_lookup_path(GrowingOnlyClustererResult& result) {
    if (graph_cache == nullptr || active_detectors.empty()) {
        return false;
    }

    auto start = Clock::now();
    std::vector<std::vector<uint32_t>> components;

    const bool buffer_covers_all_finite_pairs =
        graph_cache->has_all_pairs_distances && graph_cache->all_interior_pairs_finite &&
        graph_cache->max_finite_interior_distance < GROWING_ONLY_INF_DISTANCE &&
        bounds.buffer_bound >= graph_cache->max_finite_interior_distance;
    const bool force_component_id_grouping = force_at_max_level() &&
        graph_cache->interior_component_id_by_vertex.size() == graph_cache->num_nodes;
    const bool buffer_covers_each_interior_component =
        ((graph_cache->interior_diameter_upper_bound < GROWING_ONLY_INF_DISTANCE &&
          bounds.buffer_bound >= graph_cache->interior_diameter_upper_bound) || force_component_id_grouping) &&
        graph_cache->interior_component_id_by_vertex.size() == graph_cache->num_nodes;

    if (buffer_covers_all_finite_pairs) {
        components.emplace_back();
        components.back().reserve(active_detectors.size());
        for (size_t i = 0; i < active_detectors.size(); i++) {
            components.back().push_back(static_cast<uint32_t>(i));
        }
    } else if (buffer_covers_each_interior_component) {
        std::unordered_map<size_t, size_t> component_id_to_output_index;
        components.reserve(active_detectors.size());
        for (size_t i = 0; i < active_detectors.size(); i++) {
            auto detector = static_cast<size_t>(active_detectors[i]);
            size_t component_id = detector < graph_cache->interior_component_id_by_vertex.size()
                                      ? graph_cache->interior_component_id_by_vertex[detector]
                                      : SIZE_MAX;
            if (component_id == SIZE_MAX) {
                components.push_back({static_cast<uint32_t>(i)});
                continue;
            }
            auto it = component_id_to_output_index.find(component_id);
            if (it == component_id_to_output_index.end()) {
                size_t out = components.size();
                component_id_to_output_index.emplace(component_id, out);
                components.emplace_back();
                components.back().push_back(static_cast<uint32_t>(i));
            } else {
                components[it->second].push_back(static_cast<uint32_t>(i));
            }
        }
    } else if (graph_cache->has_radius_neighbor_precompute &&
               level < graph_cache->radius_neighbors_by_level.size()) {
        const auto& radius_neighbors = graph_cache->radius_neighbors_by_level[level];
        if (radius_neighbors.buffer_bound != bounds.buffer_bound || radius_neighbors.offsets.empty()) {
            return false;
        }

        static thread_local std::vector<uint32_t> detector_active_stamp;
        static thread_local std::vector<uint32_t> detector_active_index;
        static thread_local uint32_t detector_active_epoch = 1;
        static thread_local std::vector<uint8_t> active_visited;
        static thread_local std::vector<uint32_t> stack;

        if (detector_active_stamp.size() != graph_cache->num_nodes) {
            detector_active_stamp.assign(graph_cache->num_nodes, 0);
            detector_active_index.resize(graph_cache->num_nodes);
            detector_active_epoch = 1;
        } else if (detector_active_epoch == std::numeric_limits<uint32_t>::max()) {
            std::fill(detector_active_stamp.begin(), detector_active_stamp.end(), 0);
            detector_active_epoch = 1;
        } else {
            detector_active_epoch++;
        }
        for (size_t k = 0; k < active_detectors.size(); k++) {
            auto detector = static_cast<size_t>(active_detectors[k]);
            if (detector < detector_active_stamp.size()) {
                detector_active_stamp[detector] = detector_active_epoch;
                detector_active_index[detector] = static_cast<uint32_t>(k);
            }
        }

        active_visited.assign(active_detectors.size(), 0);
        stack.clear();
        components.reserve(std::min<size_t>(active_detectors.size(), 32));
        for (size_t start_index = 0; start_index < active_detectors.size(); start_index++) {
            if (active_visited[start_index]) {
                continue;
            }
            components.emplace_back();
            auto& component = components.back();
            active_visited[start_index] = 1;
            stack.push_back(static_cast<uint32_t>(start_index));
            while (!stack.empty()) {
                uint32_t i = stack.back();
                stack.pop_back();
                component.push_back(i);
                auto detector = static_cast<size_t>(active_detectors[i]);
                if (detector + 1 >= radius_neighbors.offsets.size()) {
                    continue;
                }
                uint64_t begin = radius_neighbors.offsets[detector];
                uint64_t end = radius_neighbors.offsets[detector + 1];
                for (uint64_t p = begin; p < end; p++) {
                    uint32_t neighbor = radius_neighbors.neighbors[p];
                    if (neighbor >= detector_active_stamp.size() ||
                        detector_active_stamp[neighbor] != detector_active_epoch) {
                        continue;
                    }
                    uint32_t j = detector_active_index[neighbor];
                    if (active_visited[j]) {
                        continue;
                    }
                    active_visited[j] = 1;
                    stack.push_back(j);
                }
            }
            std::sort(component.begin(), component.end());
        }
    } else {
        return false;
    }
    stats.event_generation_wall_ns += elapsed_ns(start);
    stats.used_component_lookup = true;

    auto processing_start = Clock::now();
    result.accepted_clusters.clear();
    result.residual_active_detectors.clear();
    result.accepted_clusters.reserve(components.size());
    const bool forced = force_at_max_level();
    const bool skip_forced_exact = skip_forced_max_exact_metadata();
    const size_t max_active_detectors = max_active_detectors_for_this_level();

    auto nearest_boundary_for_component = [&](const std::vector<uint32_t>& component) {
        cumulative_time_int best = GROWING_ONLY_INF_DISTANCE;
        for (uint32_t active_index : component) {
            best = std::min(best, nearest_boundary_match_distance_for_detector(active_detectors[active_index]));
        }
        return best;
    };

    auto farthest_boundary_for_component = [&](const std::vector<uint32_t>& component) {
        cumulative_time_int worst = 0;
        for (uint32_t active_index : component) {
            worst = std::max(worst, nearest_boundary_match_distance_for_detector(active_detectors[active_index]));
        }
        return worst;
    };

    auto exact_diameter_for_component = [&](const std::vector<uint32_t>& component, bool* exceeded_bound) {
        cumulative_time_int diameter = 0;
        *exceeded_bound = false;
        if (component.size() < 2) {
            return diameter;
        }
        stats.exact_diameter_checks++;
        for (size_t a = 0; a < component.size() && !*exceeded_bound; a++) {
            std::vector<cumulative_time_int> distances;
            if (graph_cache == nullptr || !graph_cache->has_all_pairs_distances) {
                distances = dijkstra_from_source(active_detectors[component[a]]);
            }
            for (size_t b = a + 1; b < component.size(); b++) {
                cumulative_time_int distance;
                if (graph_cache != nullptr && graph_cache->has_all_pairs_distances) {
                    distance = distance_between_active_detectors(component[a], component[b]);
                } else {
                    distance = distances[active_detectors[component[b]]];
                }
                stats.exact_internal_pair_checks++;
                if (distance >= GROWING_ONLY_INF_DISTANCE) {
                    diameter = GROWING_ONLY_INF_DISTANCE;
                    *exceeded_bound = !forced;
                    break;
                }
                diameter = std::max(diameter, distance);
                stats.internal_collisions_recorded++;
                if (!forced && distance > bounds.diameter_bound) {
                    *exceeded_bound = true;
                    break;
                }
            }
        }
        return diameter;
    };

    auto emit_component = [&](const std::vector<uint32_t>& component,
                              cumulative_time_int diameter,
                              cumulative_time_int nearest_boundary,
                              bool forced_by_max_level) {
        GrowingOnlyAcceptedCluster cluster;
        cluster.active_detectors.reserve(component.size());
        for (uint32_t active_index : component) {
            cluster.active_detectors.push_back(active_detectors[active_index]);
        }
        // component indices are sorted in the radius-neighbor path and follow
        // active-detector order in the all-finite/component-id paths.
        cluster.ready_time = 0;
        cluster.max_internal_pair_collision_time = 0;
        cluster.max_internal_pair_distance = diameter;
        cluster.first_external_collision_time = GROWING_ONLY_INF_DISTANCE;
        cluster.nearest_boundary_match_distance = nearest_boundary;
        cluster.parity_can_stop_locally = true;
        cluster.forced_by_max_level = forced_by_max_level;
        cluster.diameter = diameter;
        cluster.buffer = GROWING_ONLY_INF_DISTANCE;
        stats.max_component_active_size = std::max(stats.max_component_active_size, cluster.active_detectors.size());
        stats.removed_or_hidden_clusters++;
        if (forced_by_max_level) {
            stats.forced_by_max_level_clusters++;
            stats.forced_max_level_shortcut_count++;
        }
        result.accepted_clusters.push_back(std::move(cluster));
    };

    auto add_to_residual = [&](const std::vector<uint32_t>& component) {
        for (uint32_t active_index : component) {
            result.residual_active_detectors.push_back(active_detectors[active_index]);
        }
    };

    for (const auto& component : components) {
        if (component.empty()) {
            continue;
        }
        stats.max_component_active_size = std::max(stats.max_component_active_size, component.size());
        if (forced && skip_forced_exact) {
            emit_component(component, 0, GROWING_ONLY_INF_DISTANCE, true);
            continue;
        }
        if (!forced && max_active_detectors != 0 && component.size() > max_active_detectors) {
            stats.max_active_rejected_clusters++;
            add_to_residual(component);
            continue;
        }

        cumulative_time_int nearest_boundary = GROWING_ONLY_INF_DISTANCE;
        bool parity_ok = forced || ((component.size() & 1) == 0);
        if (!parity_ok) {
            nearest_boundary = nearest_boundary_for_component(component);
            const auto boundary_inclusive_diameter = farthest_boundary_for_component(component);
            parity_ok = boundary_inclusive_diameter <= bounds.diameter_bound;
        }
        if (!parity_ok) {
            stats.parity_rejected_clusters++;
            add_to_residual(component);
            continue;
        }

        bool exceeded_bound = false;
        cumulative_time_int diameter = exact_diameter_for_component(component, &exceeded_bound);
        if (!forced && exceeded_bound) {
            stats.diameter_rejected_clusters++;
            add_to_residual(component);
            continue;
        }
        emit_component(component, diameter, nearest_boundary, false);
    }

    if (!result.residual_active_detectors.empty()) {
        std::sort(result.residual_active_detectors.begin(), result.residual_active_detectors.end());
        result.residual_active_detectors.erase(
            std::unique(result.residual_active_detectors.begin(), result.residual_active_detectors.end()),
            result.residual_active_detectors.end());
    }
    stats.processed_pair_collision_events = stats.exact_internal_pair_checks + stats.exact_external_pair_checks;
    stats.event_processing_wall_ns += elapsed_ns(processing_start);
    result.stats = stats;
    return true;
}

void GrowingOnlyClusterer::run_sparse_frontier_streaming_unions() {
    if (active_detectors.size() < 2) {
        return;
    }

    struct FrontierEntry {
        cumulative_time_int distance = 0;
        uint32_t node = 0;
        uint32_t source = 0;
        bool operator>(const FrontierEntry& other) const {
            if (distance != other.distance) return distance > other.distance;
            if (node != other.node) return node > other.node;
            return source > other.source;
        }
    };
    struct CollisionGreater {
        bool operator()(const PairCollisionEvent& a, const PairCollisionEvent& b) const {
            return b < a;
        }
    };

    static thread_local std::vector<uint32_t> visited_epoch;
    static thread_local std::vector<cumulative_time_int> distances;
    static thread_local std::vector<uint32_t> owner;
    static thread_local uint32_t current_epoch = 1;

    if (visited_epoch.size() != graph.nodes.size()) {
        visited_epoch.assign(graph.nodes.size(), 0);
        distances.resize(graph.nodes.size());
        owner.resize(graph.nodes.size());
        current_epoch = 1;
    } else if (current_epoch == std::numeric_limits<uint32_t>::max()) {
        std::fill(visited_epoch.begin(), visited_epoch.end(), 0);
        current_epoch = 1;
    } else {
        current_epoch++;
    }

    std::priority_queue<FrontierEntry, std::vector<FrontierEntry>, std::greater<FrontierEntry>> frontier_queue;
    std::priority_queue<PairCollisionEvent, std::vector<PairCollisionEvent>, CollisionGreater> collision_queue;

    for (size_t source_index = 0; source_index < active_detectors.size(); source_index++) {
        auto source = static_cast<size_t>(active_detectors[source_index]);
        if (source >= graph.nodes.size()) {
            continue;
        }
        visited_epoch[source] = current_epoch;
        distances[source] = 0;
        owner[source] = static_cast<uint32_t>(source_index);
        frontier_queue.push({0, static_cast<uint32_t>(source), static_cast<uint32_t>(source_index)});
    }

    auto process_ready_collisions = [&](cumulative_time_int frontier_distance_floor) {
        while (!collision_queue.empty() && collision_queue.top().distance <= frontier_distance_floor) {
            auto event = collision_queue.top();
            collision_queue.pop();
            stats.processed_pair_collision_events++;
            stats.processed_collision_events++;
            auto ra = find(event.a);
            auto rb = find(event.b);
            if (meta[ra].removed_or_hidden || meta[rb].removed_or_hidden) {
                continue;
            }
            if (ra == rb) {
                continue;
            }
            unite(ra, rb, event.time, event.distance);
        }
    };

    while (!frontier_queue.empty()) {
        process_ready_collisions(frontier_queue.top().distance);
        auto entry = frontier_queue.top();
        frontier_queue.pop();
        stats.processed_frontier_events++;
        if (entry.node >= visited_epoch.size() || visited_epoch[entry.node] != current_epoch ||
            owner[entry.node] != entry.source || distances[entry.node] != entry.distance) {
            continue;
        }
        if (entry.distance > bounds.buffer_bound) {
            continue;
        }
        const DetectorNode& detector_node = graph.nodes[entry.node];
        for (size_t k = 0; k < detector_node.neighbors.size(); k++) {
            const DetectorNode* neighbor = detector_node.neighbors[k];
            if (neighbor == nullptr) {
                continue;
            }
            auto neighbor_index = static_cast<size_t>(neighbor - graph.nodes.data());
            auto edge_weight = static_cast<cumulative_time_int>(detector_node.neighbor_weights[k]);
            if (entry.distance > GROWING_ONLY_INF_DISTANCE - edge_weight) {
                continue;
            }

            if (neighbor_index < visited_epoch.size() && visited_epoch[neighbor_index] == current_epoch &&
                owner[neighbor_index] != entry.source) {
                if (entry.distance <= GROWING_ONLY_INF_DISTANCE - edge_weight &&
                    entry.distance + edge_weight <= GROWING_ONLY_INF_DISTANCE - distances[neighbor_index]) {
                    auto witness_distance = entry.distance + edge_weight + distances[neighbor_index];
                    if (witness_distance <= bounds.buffer_bound) {
                        uint32_t a = entry.source;
                        uint32_t b = owner[neighbor_index];
                        if (a > b) {
                            std::swap(a, b);
                        }
                        collision_queue.push(PairCollisionEvent{
                            static_cast<cumulative_time_int>(witness_distance / 2),
                            witness_distance,
                            a,
                            b,
                        });
                        stats.generated_pair_collision_events++;
                        stats.generated_collision_events++;
                        stats.sparse_frontier_events++;
                    }
                }
            }

            auto candidate = entry.distance + edge_weight;
            if (candidate > bounds.buffer_bound) {
                continue;
            }
            if (visited_epoch[neighbor_index] != current_epoch || candidate < distances[neighbor_index] ||
                (candidate == distances[neighbor_index] && entry.source < owner[neighbor_index])) {
                visited_epoch[neighbor_index] = current_epoch;
                distances[neighbor_index] = candidate;
                owner[neighbor_index] = entry.source;
                frontier_queue.push({candidate, static_cast<uint32_t>(neighbor_index), entry.source});
            }
        }
    }
    process_ready_collisions(GROWING_ONLY_INF_DISTANCE);
}

bool GrowingOnlyClusterer::run_sparse_frontier_path(GrowingOnlyClustererResult& result) {
    stats.used_sparse_frontier = true;
    auto generation_start = Clock::now();
    run_sparse_frontier_streaming_unions();
    stats.event_generation_wall_ns += elapsed_ns(generation_start);

    auto processing_start = Clock::now();
    // Existing preprocessing validates clusters using exact diameter when needed.
    // Keep correctness by doing this pass only for final acceptance candidates;
    // forced max-level components deliberately skip exact all-pairs metadata.
    compute_exact_component_pair_metadata();
    if (compute_exact_external_metadata()) {
        compute_first_external_collision_times();
    }
    evaluate_final_components(result);
    stats.event_processing_wall_ns += elapsed_ns(processing_start);
    result.stats = stats;
    collect_residual(result);
    return true;
}

uint32_t GrowingOnlyClusterer::find(uint32_t k) {
    while (meta[k].uf_parent != k) {
        meta[k].uf_parent = meta[meta[k].uf_parent].uf_parent;
        k = meta[k].uf_parent;
    }
    return k;
}

uint32_t GrowingOnlyClusterer::unite(
    uint32_t a,
    uint32_t b,
    cumulative_time_int event_time,
    cumulative_time_int event_distance) {
    a = find(a);
    b = find(b);
    if (a == b) {
        note_internal_collision(a, event_time, event_distance);
        return a;
    }
    if (meta[a].size < meta[b].size) {
        std::swap(a, b);
    }
    meta[b].uf_parent = a;
    meta[a].size += meta[b].size;
    meta[a].ready_time = std::max(std::max(meta[a].ready_time, meta[b].ready_time), event_time);
    meta[a].max_internal_pair_collision_time = std::max(
        std::max(meta[a].max_internal_pair_collision_time, meta[b].max_internal_pair_collision_time),
        event_time);
    meta[a].max_internal_pair_distance = std::max(
        std::max(meta[a].max_internal_pair_distance, meta[b].max_internal_pair_distance),
        event_distance);
    meta[a].nearest_boundary_match_distance = std::min(
        meta[a].nearest_boundary_match_distance, meta[b].nearest_boundary_match_distance);
    // After a merge, the newly formed connected component's buffer must be
    // measured from this new ready_time, so earlier external contacts of its
    // children are deliberately discarded.
    meta[a].first_external_collision_time = GROWING_ONLY_INF_DISTANCE;
    meta[a].accepted = false;
    meta[a].removed_or_hidden = false;

    members_by_root[a].insert(members_by_root[a].end(), members_by_root[b].begin(), members_by_root[b].end());
    members_by_root[b].clear();
    internal_collision_count_by_root[a] += internal_collision_count_by_root[b];
    external_collision_count_by_root[a] += external_collision_count_by_root[b];
    stats.union_count++;
    return a;
}

void GrowingOnlyClusterer::note_internal_collision(
    uint32_t root,
    cumulative_time_int event_time,
    cumulative_time_int event_distance) {
    root = find(root);
    meta[root].max_internal_pair_collision_time = std::max(meta[root].max_internal_pair_collision_time, event_time);
    meta[root].max_internal_pair_distance = std::max(meta[root].max_internal_pair_distance, event_distance);
    internal_collision_count_by_root[root]++;
    stats.internal_collisions_recorded++;
}

void GrowingOnlyClusterer::note_external_collision(uint32_t root, cumulative_time_int event_time) {
    root = find(root);
    if (meta[root].removed_or_hidden || event_time < meta[root].ready_time) {
        return;
    }
    meta[root].first_external_collision_time =
        std::min(meta[root].first_external_collision_time, event_time);
    external_collision_count_by_root[root]++;
    stats.external_collisions_recorded++;
}

cumulative_time_int GrowingOnlyClusterer::nearest_boundary_match_distance_for_detector(uint64_t detector) const {
    if (graph_cache != nullptr && detector < graph_cache->nearest_boundary_match_distance_by_vertex.size()) {
        return graph_cache->nearest_boundary_match_distance_by_vertex[detector];
    }

    auto distances = dijkstra_from_source(detector);
    cumulative_time_int best = GROWING_ONLY_INF_DISTANCE;
    for (size_t node = 0; node < graph.nodes.size(); node++) {
        if (distances[node] >= GROWING_ONLY_INF_DISTANCE) {
            continue;
        }
        const auto& detector_node = graph.nodes[node];
        for (size_t k = 0; k < detector_node.neighbors.size(); k++) {
            if (detector_node.neighbors[k] != nullptr) {
                continue;
            }
            auto edge_weight = static_cast<cumulative_time_int>(detector_node.neighbor_weights[k]);
            if (distances[node] > GROWING_ONLY_INF_DISTANCE - edge_weight) {
                continue;
            }
            best = std::min(best, distances[node] + edge_weight);
        }
    }
    return best;
}

cumulative_time_int GrowingOnlyClusterer::farthest_boundary_match_distance_for_root(uint32_t root) {
    root = find(root);
    cumulative_time_int worst = 0;
    for (uint32_t active_index : members_by_root[root]) {
        worst = std::max(worst, nearest_boundary_match_distance_for_detector(active_detectors[active_index]));
    }
    return worst;
}

bool GrowingOnlyClusterer::parity_can_stop_locally(uint32_t root) {
    root = find(root);
    if ((meta[root].size & 1) == 0) {
        return true;
    }
    return farthest_boundary_match_distance_for_root(root) <= bounds.diameter_bound;
}

bool GrowingOnlyClusterer::force_at_max_level() const {
    return config.force_cluster_at_max_level && level >= processing_cluster_effective_max_level(config);
}

bool GrowingOnlyClusterer::forced_phi_buffer_certificate_allows_component_split() const {
    if (!force_at_max_level()) {
        return true;
    }
    return processing_cluster_bounds_phi_buffer_certificate_ok(bounds);
}

bool GrowingOnlyClusterer::run_forced_max_level_single_root_certificate_fallback(
    GrowingOnlyClustererResult& result) {
    if (!force_at_max_level() || active_detectors.empty() ||
        forced_phi_buffer_certificate_allows_component_split()) {
        return false;
    }

    auto start = Clock::now();
    result.accepted_clusters.clear();
    result.residual_active_detectors.clear();

    GrowingOnlyAcceptedCluster cluster;
    cluster.active_detectors = std::move(active_detectors);
    std::sort(cluster.active_detectors.begin(), cluster.active_detectors.end());
    cluster.ready_time = 0;
    cluster.max_internal_pair_collision_time = 0;
    cluster.max_internal_pair_distance = 0;
    cluster.first_external_collision_time = GROWING_ONLY_INF_DISTANCE;
    cluster.nearest_boundary_match_distance = GROWING_ONLY_INF_DISTANCE;
    cluster.parity_can_stop_locally = true;
    cluster.forced_by_max_level = true;
    cluster.diameter = 0;
    cluster.buffer = GROWING_ONLY_INF_DISTANCE;

    stats.max_component_active_size = std::max(stats.max_component_active_size, cluster.active_detectors.size());
    stats.max_residual_active_size = std::max(stats.max_residual_active_size, cluster.active_detectors.size());
    stats.forced_by_max_level_clusters++;
    stats.forced_max_level_shortcut_count++;
    stats.removed_or_hidden_clusters++;
    stats.event_generation_wall_ns += elapsed_ns(start);

    result.accepted_clusters.push_back(std::move(cluster));
    result.stats = stats;
    return true;
}

size_t GrowingOnlyClusterer::max_active_detectors_for_this_level() const {
    if (level == 0 || level > config.max_active_detectors_by_level.size()) {
        return 0;
    }
    return config.max_active_detectors_by_level[level - 1];
}

void GrowingOnlyClusterer::evaluate_final_components(GrowingOnlyClustererResult& result) {
    const bool forced = force_at_max_level();
    const size_t max_active_detectors = max_active_detectors_for_this_level();
    for (size_t i = 0; i < meta.size(); i++) {
        auto root = find(static_cast<uint32_t>(i));
        if (root != i || meta[root].removed_or_hidden || members_by_root[root].empty()) {
            continue;
        }
        if (!forced && max_active_detectors != 0 && meta[root].size > max_active_detectors) {
            stats.max_active_rejected_clusters++;
            continue;
        }
        const auto diameter = meta[root].max_internal_pair_distance;
        if (!(diameter <= bounds.diameter_bound || forced)) {
            stats.diameter_rejected_clusters++;
            continue;
        }
        if (!forced && !parity_can_stop_locally(root)) {
            stats.parity_rejected_clusters++;
            continue;
        }
        emit_accepted_cluster(root, result, forced && (diameter > bounds.diameter_bound || skip_forced_max_exact_metadata()));
    }
}

void GrowingOnlyClusterer::emit_accepted_cluster(
    uint32_t root,
    GrowingOnlyClustererResult& result,
    bool forced_by_max_level) {
    root = find(root);
    if (meta[root].removed_or_hidden) {
        return;
    }
    meta[root].accepted = true;
    meta[root].removed_or_hidden = true;
    GrowingOnlyAcceptedCluster cluster;
    cluster.active_detectors.reserve(members_by_root[root].size());
    for (auto active_index : members_by_root[root]) {
        cluster.active_detectors.push_back(active_detectors[active_index]);
    }
    stats.max_component_active_size = std::max(stats.max_component_active_size, cluster.active_detectors.size());
    std::sort(cluster.active_detectors.begin(), cluster.active_detectors.end());
    cluster.ready_time = meta[root].ready_time;
    cluster.max_internal_pair_collision_time = meta[root].max_internal_pair_collision_time;
    cluster.max_internal_pair_distance = meta[root].max_internal_pair_distance;
    cluster.first_external_collision_time = meta[root].first_external_collision_time;
    cluster.nearest_boundary_match_distance = meta[root].nearest_boundary_match_distance;
    cluster.parity_can_stop_locally = force_at_max_level() || parity_can_stop_locally(root);
    cluster.forced_by_max_level = forced_by_max_level;
    cluster.diameter = cluster.max_internal_pair_distance;
    cluster.buffer = cluster.first_external_collision_time >= GROWING_ONLY_INF_DISTANCE
                         ? GROWING_ONLY_INF_DISTANCE
                         : cluster.first_external_collision_time - cluster.ready_time;
    cluster.internal_collisions_recorded = internal_collision_count_by_root[root];
    cluster.external_collisions_recorded = external_collision_count_by_root[root];
    result.accepted_clusters.push_back(std::move(cluster));
    stats.removed_or_hidden_clusters++;
    if (forced_by_max_level) {
        stats.forced_by_max_level_clusters++;
    }
}

void GrowingOnlyClusterer::collect_residual(GrowingOnlyClustererResult& result) {
    result.residual_active_detectors.clear();
    for (size_t i = 0; i < active_detectors.size(); i++) {
        auto root = find(static_cast<uint32_t>(i));
        if (meta[root].removed_or_hidden) {
            continue;
        }
        result.residual_active_detectors.push_back(active_detectors[i]);
    }
    std::sort(result.residual_active_detectors.begin(), result.residual_active_detectors.end());
}


bool GrowingOnlyClusterer::try_run_no_union_singleton_fast_path(GrowingOnlyClustererResult& result) {
    if (active_detectors.empty()) {
        return false;
    }
    auto start = Clock::now();

    auto accept_as_singletons_direct = [&]() {
        stats.event_generation_wall_ns += elapsed_ns(start);
        auto processing_start = Clock::now();
        result.accepted_clusters.clear();
        result.residual_active_detectors.clear();
        const bool forced = force_at_max_level();
        const size_t max_active_detectors = max_active_detectors_for_this_level();

        // In the common match-path selected policy, level 1 proves that no
        // detector-detector union is possible and every singleton is odd.  If
        // even the graph-wide nearest boundary distance is above the diameter
        // bound, all singletons are guaranteed parity-rejected.  Preserve the
        // exact decision while avoiding one nearest-boundary lookup and one
        // push_back per active detector.
        if (!forced && max_active_detectors != 0 && max_active_detectors < 1) {
            stats.max_active_rejected_clusters += active_detectors.size();
            result.residual_active_detectors = std::move(active_detectors);
            stats.event_processing_wall_ns += elapsed_ns(processing_start);
            result.stats = stats;
            return true;
        }
        if (!forced && cached_min_nearest_boundary_match_distance(graph_cache) > bounds.diameter_bound) {
            stats.parity_rejected_clusters += active_detectors.size();
            result.residual_active_detectors = std::move(active_detectors);
            stats.event_processing_wall_ns += elapsed_ns(processing_start);
            result.stats = stats;
            return true;
        }

        result.accepted_clusters.reserve(active_detectors.size());
        for (auto detector : active_detectors) {
            if (!forced && max_active_detectors != 0 && max_active_detectors < 1) {
                result.residual_active_detectors.push_back(detector);
                stats.max_active_rejected_clusters++;
                continue;
            }
            cumulative_time_int nearest_boundary = GROWING_ONLY_INF_DISTANCE;
            bool parity_ok = forced;
            if (!forced) {
                nearest_boundary = nearest_boundary_match_distance_for_detector(detector);
                parity_ok = nearest_boundary <= bounds.diameter_bound;
            }
            if (!parity_ok) {
                result.residual_active_detectors.push_back(detector);
                stats.parity_rejected_clusters++;
                continue;
            }

            GrowingOnlyAcceptedCluster cluster;
            cluster.active_detectors.push_back(detector);
            cluster.ready_time = 0;
            cluster.max_internal_pair_collision_time = 0;
            cluster.max_internal_pair_distance = 0;
            cluster.first_external_collision_time = GROWING_ONLY_INF_DISTANCE;
            cluster.nearest_boundary_match_distance = nearest_boundary;
            cluster.parity_can_stop_locally = true;
            cluster.forced_by_max_level = forced;
            cluster.diameter = 0;
            cluster.buffer = GROWING_ONLY_INF_DISTANCE;
            result.accepted_clusters.push_back(std::move(cluster));
            stats.max_component_active_size = std::max<size_t>(stats.max_component_active_size, 1);
            stats.removed_or_hidden_clusters++;
            if (forced) {
                stats.forced_by_max_level_clusters++;
            }
        }
        stats.processed_pair_collision_events = 0;
        stats.event_processing_wall_ns += elapsed_ns(processing_start);
        result.stats = stats;
        return true;
    };

    // Cheap graph-wide certificate: no detector-detector union can happen when
    // the level buffer is smaller than every interior edge. This is safe whether
    // or not all-pairs distances are available and avoids dense active-pair
    // singleton scans and UF materialization.
    if (bounds.buffer_bound < cached_min_interior_edge_weight(graph)) {
        return accept_as_singletons_direct();
    }

    // Default sparse singleton certificate.  The radius-neighbor lookup is a
    // graph-only precompute keyed by the current buffer bound; checking if any
    // active detector has an active neighbor inside the bound is O(k + touched
    // radius-neighbors), not O(k^2).  If there is no such pair, the singleton
    // fast path is valid.  If there is one, fall through to component-lookup /
    // sparse-frontier clustering.
    if (graph_cache != nullptr && graph_cache->has_radius_neighbor_precompute &&
        level < graph_cache->radius_neighbors_by_level.size()) {
        const auto& radius_neighbors = graph_cache->radius_neighbors_by_level[level];
        if (radius_neighbors.buffer_bound == bounds.buffer_bound && !radius_neighbors.offsets.empty()) {
            static thread_local std::vector<uint32_t> detector_active_stamp;
            static thread_local uint32_t detector_active_epoch = 1;
            if (detector_active_stamp.size() != graph_cache->num_nodes) {
                detector_active_stamp.assign(graph_cache->num_nodes, 0);
                detector_active_epoch = 1;
            } else if (detector_active_epoch == std::numeric_limits<uint32_t>::max()) {
                std::fill(detector_active_stamp.begin(), detector_active_stamp.end(), 0);
                detector_active_epoch = 1;
            } else {
                detector_active_epoch++;
            }
            for (auto detector_u64 : active_detectors) {
                auto detector = static_cast<size_t>(detector_u64);
                if (detector < detector_active_stamp.size()) {
                    detector_active_stamp[detector] = detector_active_epoch;
                }
            }
            bool saw_union_eligible_pair = false;
            for (auto detector_u64 : active_detectors) {
                auto detector = static_cast<size_t>(detector_u64);
                if (detector + 1 >= radius_neighbors.offsets.size()) {
                    continue;
                }
                uint64_t begin = radius_neighbors.offsets[detector];
                uint64_t end = radius_neighbors.offsets[detector + 1];
                for (uint64_t p = begin; p < end; p++) {
                    uint32_t neighbor = radius_neighbors.neighbors[p];
                    if (neighbor == detector) {
                        continue;
                    }
                    if (neighbor < detector_active_stamp.size() &&
                        detector_active_stamp[neighbor] == detector_active_epoch) {
                        saw_union_eligible_pair = true;
                        break;
                    }
                }
                if (saw_union_eligible_pair) {
                    break;
                }
            }
            if (!saw_union_eligible_pair) {
                return accept_as_singletons_direct();
            }
            return false;
        }
    }

    // Reference/debug fallback only.  This used to be the default when the
    // all-pairs distance table existed, but it is O(k^2) in active detectors and
    // was a hidden source of preprocessing cost.
    if (!env_is_enabled("PYMATCHING_GROWING_ONLY_SINGLETON_ALL_PAIR_CHECK")) {
        return false;
    }
    if (graph_cache == nullptr || !graph_cache->has_all_pairs_distances) {
        return false;
    }

    bool saw_union_eligible_pair = false;
    size_t pair_count = 0;
    for (size_t i = 0; i < active_detectors.size(); i++) {
        for (size_t j = i + 1; j < active_detectors.size(); j++) {
            auto distance = distance_between_active_detectors(i, j);
            if (distance >= GROWING_ONLY_INF_DISTANCE) {
                continue;
            }
            pair_count++;
            if (distance <= bounds.buffer_bound) {
                saw_union_eligible_pair = true;
                break;
            }
        }
        if (saw_union_eligible_pair) {
            break;
        }
    }
    if (saw_union_eligible_pair) {
        return false;
    }
    stats.generated_pair_collision_events = pair_count;
    return accept_as_singletons_direct();
}

GrowingOnlyClustererResult GrowingOnlyClusterer::run() {
    GrowingOnlyClustererResult result;
    result.level = level;
    result.bounds = bounds;
    stats.max_residual_active_size = active_detectors.size();

    if (run_forced_max_level_single_root_certificate_fallback(result)) {
        return result;
    }

    if (run_forced_max_level_direct_component_lookup_path(result)) {
        return result;
    }

    // The default sparse/component-lookup paths can now accept/reject directly
    // from component lists.  Try them before allocating UF metadata and
    // members_by_root; that materialization is retained only for sparse-frontier
    // and active-pair reference fallbacks.
    if (try_run_no_union_singleton_fast_path(result)) {
        return result;
    }
    if (use_sparse_frontier_path() && run_precomputed_component_lookup_path(result)) {
        return result;
    }

    meta.resize(active_detectors.size());
    members_by_root.resize(active_detectors.size());
    internal_collision_count_by_root.assign(active_detectors.size(), 0);
    external_collision_count_by_root.assign(active_detectors.size(), 0);
    for (size_t i = 0; i < active_detectors.size(); i++) {
        meta[i].uf_parent = static_cast<uint32_t>(i);
        meta[i].size = 1;
        meta[i].nearest_boundary_match_distance = nearest_boundary_match_distance_for_detector(active_detectors[i]);
        members_by_root[i].push_back(static_cast<uint32_t>(i));
    }

    if (use_sparse_frontier_path()) {
        run_sparse_frontier_path(result);
        return result;
    }

    if (run_active_pair_table_component_path(result)) {
        return result;
    }

    auto generation_start = Clock::now();
    auto events = generate_pair_collision_events();
    stats.event_generation_wall_ns += elapsed_ns(generation_start);
    stats.generated_pair_collision_events = events.size();

    auto processing_start = Clock::now();
    for (const auto& event : events) {
        auto ra = find(event.a);
        auto rb = find(event.b);
        if (meta[ra].removed_or_hidden || meta[rb].removed_or_hidden) {
            continue;
        }
        stats.processed_pair_collision_events++;
        if (ra == rb) {
            // Same-component collisions after a chain has already connected the
            // component are kept. They are required to recover the existing
            // ProcessingCluster diameter check, which uses the maximum pair
            // distance inside the buffer-connected component, not just merge
            // edges.
            note_internal_collision(ra, event.time, event.distance);
            continue;
        }

        if (event.distance <= bounds.buffer_bound) {
            // Existing preprocessing constructs components using the level's
            // buffer bound as a detector-graph distance threshold. Mirror that
            // rule here instead of the previous diameter-bound union heuristic.
            unite(ra, rb, event.time, event.distance);
        } else {
            // This collision is between final buffer-components. It is useful
            // metadata for experimental buffer reporting, but acceptance is
            // decided only after all internal collisions have been observed.
            note_external_collision(ra, event.time);
            note_external_collision(rb, event.time);
        }
    }
    evaluate_final_components(result);
    stats.event_processing_wall_ns += elapsed_ns(processing_start);
    result.stats = stats;
    collect_residual(result);
    return result;
}



cumulative_time_int cached_all_pairs_distance_for_adaptive_bounds(
    const ProcessingClusterGraphCache* graph_cache,
    uint64_t source,
    uint64_t destination) {
    if (graph_cache == nullptr || !graph_cache->has_all_pairs_distances ||
        source >= graph_cache->num_nodes || destination >= graph_cache->num_nodes) {
        return GROWING_ONLY_INF_DISTANCE;
    }
    size_t a = static_cast<size_t>(source);
    size_t b = static_cast<size_t>(destination);
    if (a > b) {
        std::swap(a, b);
    }
    if (a >= graph_cache->packed_all_pairs_row_offsets.size()) {
        return GROWING_ONLY_INF_DISTANCE;
    }
    uint64_t index = graph_cache->packed_all_pairs_row_offsets[a] + (b - a);
    if (index >= graph_cache->packed_all_pairs_interior_distances.size()) {
        return GROWING_ONLY_INF_DISTANCE;
    }
    uint32_t packed = graph_cache->packed_all_pairs_interior_distances[index];
    if (packed == std::numeric_limits<uint32_t>::max()) {
        return GROWING_ONLY_INF_DISTANCE;
    }
    return static_cast<cumulative_time_int>(packed);
}


cumulative_time_int shortest_detector_distance_for_level1_short_pair(
    const MatchingGraph& graph,
    const ProcessingClusterGraphCache* graph_cache,
    uint64_t source,
    uint64_t destination,
    cumulative_time_int stop_after = GROWING_ONLY_INF_DISTANCE) {
    if (source == destination) {
        return 0;
    }
    cumulative_time_int cached = cached_all_pairs_distance_for_adaptive_bounds(graph_cache, source, destination);
    if (cached < GROWING_ONLY_INF_DISTANCE) {
        return cached;
    }
    if (source >= graph.nodes.size() || destination >= graph.nodes.size()) {
        return GROWING_ONLY_INF_DISTANCE;
    }
    std::vector<cumulative_time_int> distances(graph.nodes.size(), GROWING_ONLY_INF_DISTANCE);
    using QueueEntry = std::pair<cumulative_time_int, size_t>;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
    distances[source] = 0;
    queue.push({0, static_cast<size_t>(source)});
    while (!queue.empty()) {
        auto [distance, node] = queue.top();
        queue.pop();
        if (distance != distances[node]) {
            continue;
        }
        if (node == static_cast<size_t>(destination)) {
            return distance;
        }
        if (distance > stop_after) {
            continue;
        }
        const auto& detector_node = graph.nodes[node];
        for (size_t k = 0; k < detector_node.neighbors.size(); k++) {
            const DetectorNode* neighbor = detector_node.neighbors[k];
            if (neighbor == nullptr) {
                continue;
            }
            auto neighbor_index = static_cast<size_t>(neighbor - graph.nodes.data());
            auto edge_weight = static_cast<cumulative_time_int>(detector_node.neighbor_weights[k]);
            if (distance > GROWING_ONLY_INF_DISTANCE - edge_weight) {
                continue;
            }
            auto candidate = distance + edge_weight;
            if (candidate > stop_after && candidate >= distances[destination]) {
                continue;
            }
            if (candidate < distances[neighbor_index]) {
                distances[neighbor_index] = candidate;
                queue.push({candidate, neighbor_index});
            }
        }
    }
    return distances[destination];
}

cumulative_time_int saturating_add_for_adaptive_bounds(
    cumulative_time_int a,
    cumulative_time_int b) {
    if (a >= GROWING_ONLY_INF_DISTANCE || b >= GROWING_ONLY_INF_DISTANCE ||
        a > GROWING_ONLY_INF_DISTANCE - b) {
        return GROWING_ONLY_INF_DISTANCE;
    }
    return a + b;
}

cumulative_time_int saturating_mul_for_adaptive_bounds(
    cumulative_time_int a,
    cumulative_time_int b) {
    if (a == 0 || b == 0) {
        return 0;
    }
    if (a >= GROWING_ONLY_INF_DISTANCE || b >= GROWING_ONLY_INF_DISTANCE ||
        a > GROWING_ONLY_INF_DISTANCE / b) {
        return GROWING_ONLY_INF_DISTANCE;
    }
    return a * b;
}

void add_positive_candidate_distance(
    std::vector<cumulative_time_int>& values,
    cumulative_time_int value) {
    if (value > 0 && value < GROWING_ONLY_INF_DISTANCE) {
        values.push_back(value);
    }
}

std::vector<cumulative_time_int> unique_limited_adaptive_values(
    std::vector<cumulative_time_int> values,
    size_t limit) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    if (values.size() <= limit || limit == 0) {
        return values;
    }
    std::vector<cumulative_time_int> limited;
    limited.reserve(limit);
    for (size_t k = 0; k < limit; k++) {
        size_t index = (k * (values.size() - 1)) / (limit - 1);
        limited.push_back(values[index]);
    }
    limited.erase(std::unique(limited.begin(), limited.end()), limited.end());
    return limited;
}

void ensure_adaptive_explicit_bounds_vector(ProcessingClusterConfig& config) {
    const size_t target_levels = std::max<size_t>(config.max_level, 1);
    if (config.explicit_bounds_by_level.size() >= target_levels) {
        return;
    }
    ProcessingClusterConfig base = config;
    base.adaptive_cluster_bounds = false;
    std::vector<ProcessingClusterBounds> existing = config.explicit_bounds_by_level;
    config.explicit_bounds_by_level.resize(target_levels);
    for (size_t level = 1; level <= target_levels; level++) {
        if (level <= existing.size()) {
            config.explicit_bounds_by_level[level - 1] = existing[level - 1];
        } else {
            config.explicit_bounds_by_level[level - 1] = processing_cluster_bounds_for_level(level, base);
        }
    }
}

long double adaptive_certificate_bound_value(const ProcessingClusterBounds& bounds, bool use_diameter) {
    if (bounds.has_long_double_schedule_bounds) {
        return use_diameter ? bounds.schedule_diameter_bound : bounds.schedule_buffer_bound;
    }
    return static_cast<long double>(use_diameter ? bounds.diameter_bound : bounds.buffer_bound);
}

long double adaptive_phi_for_level(const ProcessingClusterConfig& config, size_t level) {
    if (level <= 1) {
        return 1.0L;
    }
    std::vector<long double> phi(level + 1, 0.0L);
    phi[1] = 1.0L;
    for (size_t upper = 2; upper <= level; upper++) {
        long double occupied = 0.0L;
        for (size_t lower = 1; lower < upper; lower++) {
            if (lower > config.explicit_bounds_by_level.size()) {
                return -std::numeric_limits<long double>::infinity();
            }
            const auto& lower_bounds = config.explicit_bounds_by_level[lower - 1];
            const long double d_plus_1 = adaptive_certificate_bound_value(lower_bounds, true) + 1.0L;
            const long double b_minus_1 = adaptive_certificate_bound_value(lower_bounds, false) - 1.0L;
            const long double denominator =
                (phi[lower] + 1.0L) * d_plus_1 + phi[lower] * b_minus_1;
            if (!(denominator > 0.0L) || !std::isfinite(static_cast<double>(denominator))) {
                return -std::numeric_limits<long double>::infinity();
            }
            occupied += ((phi[lower] + 2.0L) * d_plus_1) / denominator;
        }
        phi[upper] = 1.0L - occupied;
    }
    return phi[level];
}

bool adaptive_raise_buffer_to_optional_lower_bounds(
    const ProcessingClusterConfig& config,
    size_t level,
    ProcessingClusterBounds& candidate) {
    if (level > 1 && config.adaptive_enforce_monotone_bounds) {
        if (level - 2 >= config.explicit_bounds_by_level.size()) {
            return false;
        }
        candidate.buffer_bound = std::max(candidate.buffer_bound,
            config.explicit_bounds_by_level[level - 2].buffer_bound);
    }
    if (config.adaptive_min_buffer_to_diameter_ratio > 0.0 &&
        (!config.adaptive_min_buffer_ratio_first_level_only || level == 1)) {
        long double required = static_cast<long double>(candidate.diameter_bound) *
            static_cast<long double>(config.adaptive_min_buffer_to_diameter_ratio);
        if (required >= static_cast<long double>(GROWING_ONLY_INF_DISTANCE)) {
            return false;
        }
        candidate.buffer_bound = std::max<cumulative_time_int>(
            candidate.buffer_bound,
            static_cast<cumulative_time_int>(std::ceil(required)));
    }
    if (config.adaptive_enforce_phi_buffer_certificate) {
        long double phi = adaptive_phi_for_level(config, level);
        if (!(phi > 0.0L) || !std::isfinite(static_cast<double>(phi))) {
            return false;
        }
        long double required = 1.0L +
            2.0L * (static_cast<long double>(candidate.diameter_bound) + 1.0L) / phi;
        if (required >= static_cast<long double>(GROWING_ONLY_INF_DISTANCE)) {
            return false;
        }
        candidate.buffer_bound = std::max<cumulative_time_int>(
            candidate.buffer_bound,
            static_cast<cumulative_time_int>(std::ceil(required)));
    }
    return true;
}


cumulative_time_int adaptive_legacy_gap_lower_bound(
    const ProcessingClusterConfig& config,
    size_t level) {
    if (level <= 1 || !config.adaptive_enforce_legacy_gap_bound) {
        return 1;
    }
    if (level - 2 >= config.explicit_bounds_by_level.size()) {
        return GROWING_ONLY_INF_DISTANCE;
    }
    const auto& previous = config.explicit_bounds_by_level[level - 2];
    return saturating_add_for_adaptive_bounds(
        saturating_mul_for_adaptive_bounds(3, previous.diameter_bound),
        saturating_mul_for_adaptive_bounds(4, previous.buffer_bound));
}

cumulative_time_int ceil_to_cumulative_time_int_for_adaptive_bounds(long double value) {
    if (!(value < static_cast<long double>(GROWING_ONLY_INF_DISTANCE)) ||
        !std::isfinite(static_cast<double>(value))) {
        return GROWING_ONLY_INF_DISTANCE;
    }
    if (value <= 0.0L) {
        return 0;
    }
    return static_cast<cumulative_time_int>(std::ceil(value));
}

cumulative_time_int adaptive_required_buffer_for_diameter(
    const ProcessingClusterConfig& config,
    size_t level,
    cumulative_time_int diameter) {
    cumulative_time_int required = 2;
    if (level > 1 && config.adaptive_enforce_monotone_bounds) {
        if (level - 2 >= config.explicit_bounds_by_level.size()) {
            return GROWING_ONLY_INF_DISTANCE;
        }
        required = std::max(required, config.explicit_bounds_by_level[level - 2].buffer_bound);
    }
    if (config.adaptive_min_buffer_to_diameter_ratio > 0.0 &&
        (!config.adaptive_min_buffer_ratio_first_level_only || level == 1)) {
        long double ratio_required = static_cast<long double>(diameter) *
            static_cast<long double>(config.adaptive_min_buffer_to_diameter_ratio);
        required = std::max(required, ceil_to_cumulative_time_int_for_adaptive_bounds(ratio_required));
    }
    if (config.adaptive_enforce_phi_buffer_certificate) {
        long double phi = adaptive_phi_for_level(config, level);
        if (!(phi > 0.0L) || !std::isfinite(static_cast<double>(phi))) {
            return GROWING_ONLY_INF_DISTANCE;
        }
        if (config.adaptive_phi_floor > 0.0 &&
            phi + 1e-18L < static_cast<long double>(config.adaptive_phi_floor)) {
            return GROWING_ONLY_INF_DISTANCE;
        }
        long double phi_required = 1.0L +
            2.0L * (static_cast<long double>(diameter) + 1.0L) / phi;
        required = std::max(required, ceil_to_cumulative_time_int_for_adaptive_bounds(phi_required));
    }
    return required;
}

cumulative_time_int adaptive_max_diameter_for_buffer(
    const ProcessingClusterConfig& config,
    size_t level,
    cumulative_time_int buffer) {
    if (!(buffer > 1) || buffer >= GROWING_ONLY_INF_DISTANCE) {
        return 0;
    }
    long double upper = static_cast<long double>(GROWING_ONLY_INF_DISTANCE - 1);
    if (config.adaptive_enforce_phi_buffer_certificate) {
        long double phi = adaptive_phi_for_level(config, level);
        if (!(phi > 0.0L) || !std::isfinite(static_cast<double>(phi))) {
            return 0;
        }
        if (config.adaptive_phi_floor > 0.0 &&
            phi + 1e-18L < static_cast<long double>(config.adaptive_phi_floor)) {
            return 0;
        }
        upper = std::min(upper, 0.5L * phi * (static_cast<long double>(buffer) - 1.0L) - 1.0L);
    }
    if (config.adaptive_min_buffer_to_diameter_ratio > 0.0 &&
        (!config.adaptive_min_buffer_ratio_first_level_only || level == 1)) {
        upper = std::min(upper,
            static_cast<long double>(buffer) /
                static_cast<long double>(config.adaptive_min_buffer_to_diameter_ratio));
    }
    if (!(upper > 0.0L) || !std::isfinite(static_cast<double>(upper))) {
        return 0;
    }
    if (upper >= static_cast<long double>(GROWING_ONLY_INF_DISTANCE)) {
        return GROWING_ONLY_INF_DISTANCE - 1;
    }
    return static_cast<cumulative_time_int>(std::floor(upper + 1e-12L));
}

void adaptive_fill_certificate_metadata(
    const ProcessingClusterConfig& config,
    size_t level,
    ProcessingClusterBounds& bounds) {
    if (config.adaptive_enforce_phi_buffer_certificate) {
        long double phi = adaptive_phi_for_level(config, level);
        if (phi > 0.0L && std::isfinite(static_cast<double>(phi))) {
            bounds.has_phi_buffer_certificate = true;
            bounds.phi = static_cast<double>(phi);
            bounds.phi_required_buffer_bound = ceil_to_cumulative_time_int_for_adaptive_bounds(
                1.0L + 2.0L * (static_cast<long double>(bounds.diameter_bound) + 1.0L) / phi);
        }
    }
    if (config.adaptive_min_buffer_to_diameter_ratio > 0.0 &&
        (!config.adaptive_min_buffer_ratio_first_level_only || level == 1)) {
        bounds.ratio_required_buffer_bound = ceil_to_cumulative_time_int_for_adaptive_bounds(
            static_cast<long double>(bounds.diameter_bound) *
            static_cast<long double>(config.adaptive_min_buffer_to_diameter_ratio));
    }
    bounds.min_required_buffer_bound = std::max<cumulative_time_int>(
        2,
        std::max(bounds.phi_required_buffer_bound,
                 std::max(bounds.phi_budget_required_buffer_bound,
                          bounds.ratio_required_buffer_bound)));
}

bool adaptive_bounds_satisfy_optional_constraints(
    const ProcessingClusterConfig& config,
    size_t level,
    const ProcessingClusterBounds& candidate) {
    if (!(candidate.diameter_bound > 0) || !(candidate.buffer_bound > 1)) {
        return false;
    }
    if (level > 1 && (config.adaptive_enforce_monotone_bounds || config.adaptive_enforce_legacy_gap_bound)) {
        if (level - 2 >= config.explicit_bounds_by_level.size()) {
            return false;
        }
        const auto& previous = config.explicit_bounds_by_level[level - 2];
        if (config.adaptive_enforce_monotone_bounds) {
            if (candidate.diameter_bound < previous.diameter_bound ||
                candidate.buffer_bound < previous.buffer_bound) {
                return false;
            }
        }
        if (config.adaptive_enforce_legacy_gap_bound) {
            cumulative_time_int required = saturating_add_for_adaptive_bounds(
                saturating_mul_for_adaptive_bounds(3, previous.diameter_bound),
                saturating_mul_for_adaptive_bounds(4, previous.buffer_bound));
            if (candidate.diameter_bound < required) {
                return false;
            }
        }
    }
    if (config.adaptive_min_buffer_to_diameter_ratio > 0.0 &&
        (!config.adaptive_min_buffer_ratio_first_level_only || level == 1)) {
        long double required = static_cast<long double>(candidate.diameter_bound) *
            static_cast<long double>(config.adaptive_min_buffer_to_diameter_ratio);
        if (static_cast<long double>(candidate.buffer_bound) + 1e-9L < required) {
            return false;
        }
    }
    if (config.adaptive_enforce_phi_buffer_certificate) {
        long double phi = adaptive_phi_for_level(config, level);
        if (!(phi > 0.0L) || !std::isfinite(static_cast<double>(phi))) {
            return false;
        }
        if (config.adaptive_phi_floor > 0.0 &&
            phi + 1e-18L < static_cast<long double>(config.adaptive_phi_floor)) {
            return false;
        }
        long double required = 1.0L +
            2.0L * (static_cast<long double>(candidate.diameter_bound) + 1.0L) / phi;
        if (static_cast<long double>(candidate.buffer_bound) + 1e-9L < required) {
            return false;
        }
    }
    return true;
}

ProcessingClusterConfig adaptive_trial_config_with_bounds(
    const ProcessingClusterConfig& config,
    size_t level,
    const ProcessingClusterBounds& bounds) {
    ProcessingClusterConfig trial = config;
    trial.adaptive_cluster_bounds = false;
    ensure_adaptive_explicit_bounds_vector(trial);
    if (level > trial.explicit_bounds_by_level.size()) {
        trial.explicit_bounds_by_level.resize(level, bounds);
        trial.max_level = std::max(trial.max_level, level);
    }
    auto populated_bounds = bounds;
    adaptive_fill_certificate_metadata(trial, level, populated_bounds);
    trial.explicit_bounds_by_level[level - 1] = populated_bounds;
    return trial;
}

std::vector<cumulative_time_int> adaptive_candidate_values_from_residual(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& active_detectors,
    const ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache,
    size_t level,
    bool for_buffer) {
    std::vector<cumulative_time_int> values;
    if (level > 0 && level <= config.explicit_bounds_by_level.size()) {
        const auto& reference = config.explicit_bounds_by_level[level - 1];
        add_positive_candidate_distance(values, reference.diameter_bound);
        add_positive_candidate_distance(values, reference.buffer_bound);
    }
    if (level > 1 && level - 2 < config.explicit_bounds_by_level.size()) {
        const auto& previous = config.explicit_bounds_by_level[level - 2];
        add_positive_candidate_distance(values, previous.diameter_bound);
        add_positive_candidate_distance(values, previous.buffer_bound);
        add_positive_candidate_distance(values, previous.diameter_bound + 1);
        add_positive_candidate_distance(values, previous.buffer_bound + 1);
        if (config.adaptive_enforce_legacy_gap_bound) {
            add_positive_candidate_distance(values, saturating_add_for_adaptive_bounds(
                saturating_mul_for_adaptive_bounds(3, previous.diameter_bound),
                saturating_mul_for_adaptive_bounds(4, previous.buffer_bound)));
        }
    }
    if (graph_cache != nullptr && !graph_cache->nearest_boundary_match_distance_by_vertex.empty()) {
        for (auto detector : active_detectors) {
            if (detector < graph_cache->nearest_boundary_match_distance_by_vertex.size()) {
                add_positive_candidate_distance(values, graph_cache->nearest_boundary_match_distance_by_vertex[detector]);
            }
        }
    } else if (!for_buffer) {
        // Boundary distances are most important for odd-component diameter tests.
        // Avoid a Dijkstra-per-detector fallback here; the fixed reference bound
        // remains available when the cheap cache was intentionally not built.
        (void)graph;
    }

    const size_t active_count = active_detectors.size();
    const size_t sample_cap = std::max<size_t>(64, config.adaptive_candidate_limit * 16);
    size_t sampled_pairs = 0;
    for (size_t i = 0; i < active_count; i++) {
        for (size_t j = i + 1; j < active_count; j++) {
            cumulative_time_int distance = cached_all_pairs_distance_for_adaptive_bounds(
                graph_cache, active_detectors[i], active_detectors[j]);
            add_positive_candidate_distance(values, distance);
            if (++sampled_pairs >= sample_cap && active_count > 128) {
                i = active_count;
                break;
            }
        }
    }

    if (values.empty()) {
        add_positive_candidate_distance(values, 1);
        if (level > 0 && level <= config.explicit_bounds_by_level.size()) {
            add_positive_candidate_distance(values, for_buffer
                ? config.explicit_bounds_by_level[level - 1].buffer_bound
                : config.explicit_bounds_by_level[level - 1].diameter_bound);
        }
    }
    const size_t axis_limit = std::max<size_t>(4,
        static_cast<size_t>(std::ceil(std::sqrt(static_cast<double>(std::max<size_t>(1, config.adaptive_candidate_limit))))));
    return unique_limited_adaptive_values(std::move(values), axis_limit);
}

size_t accepted_active_count_for_adaptive_result(const GrowingOnlyClustererResult& result) {
    size_t total = 0;
    for (const auto& cluster : result.accepted_clusters) {
        total += cluster.active_detectors.size();
    }
    return total;
}

size_t max_accepted_cluster_active_for_adaptive_result(const GrowingOnlyClustererResult& result) {
    size_t max_active = 0;
    for (const auto& cluster : result.accepted_clusters) {
        max_active = std::max(max_active, cluster.active_detectors.size());
    }
    return max_active;
}

uint64_t adaptive_event_count_proxy_for_detectors(
    const ProcessingClusterConfig& config,
    const std::vector<uint64_t>& active_detectors) {
    if (active_detectors.empty()) {
        return 0;
    }
    if (config.adaptive_event_count_callback) {
        return config.adaptive_event_count_callback(active_detectors);
    }
    // Fallback used outside event-count-only benchmarks.  It preserves ordering
    // toward smaller clusters without pretending to be an exact MWPM event count.
    return static_cast<uint64_t>((active_detectors.size() + 1) / 2);
}

AdaptiveEventTiming adaptive_event_timing_for_detectors(
    const ProcessingClusterConfig& config,
    const std::vector<uint64_t>& active_detectors,
    cumulative_time_int cutoff_algorithmic_time) {
    AdaptiveEventTiming timing;
    if (active_detectors.empty()) {
        return timing;
    }
    if (config.adaptive_event_timing_callback) {
        return config.adaptive_event_timing_callback(active_detectors, cutoff_algorithmic_time);
    }
    const uint64_t raw_events = adaptive_event_count_proxy_for_detectors(config, active_detectors);
    timing.raw_event_count = raw_events;
    timing.events_after_cutoff = raw_events;
    timing.stop_algorithmic_time = 0;
    return timing;
}

uint64_t max_accepted_cluster_events_for_adaptive_result(
    const GrowingOnlyClustererResult& result,
    const ProcessingClusterConfig& config) {
    uint64_t max_events = 0;
    for (const auto& cluster : result.accepted_clusters) {
        max_events = std::max(max_events,
            adaptive_event_count_proxy_for_detectors(config, cluster.active_detectors));
    }
    return max_events;
}

const char* adaptive_objective_label(const ProcessingClusterConfig& config) {
    if (config.adaptive_objective_minimize_ideal_events) {
        return "ideal_events";
    }
    if (config.adaptive_objective_minimize_max_events) {
        return "max_events";
    }
    if (config.adaptive_objective_minimize_max_active) {
        return "max_active";
    }
    return "cluster_count";
}

bool adaptive_cluster_should_be_deferred(
    const GrowingOnlyAcceptedCluster& cluster,
    const ProcessingClusterConfig& config) {
    if (!config.adaptive_defer_heavy_clusters || cluster.forced_by_max_level) {
        return false;
    }
    if (config.adaptive_defer_max_active != 0 &&
        cluster.active_detectors.size() > config.adaptive_defer_max_active) {
        return true;
    }
    if (config.adaptive_defer_max_events != 0) {
        const uint64_t events = adaptive_event_count_proxy_for_detectors(config, cluster.active_detectors);
        if (events > config.adaptive_defer_max_events) {
            return true;
        }
    }
    return false;
}


void apply_level1_short_pair_acceptance(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& original_active_detectors,
    GrowingOnlyClustererResult& result,
    const ProcessingClusterBounds& bounds,
    const ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache) {
    cumulative_time_int configured_max_distance = config.level1_short_pair_max_distance;
    if (configured_max_distance == 0 && config.level1_short_pair_distance_ratio > 0.0) {
        long double scaled = static_cast<long double>(bounds.diameter_bound) *
            static_cast<long double>(config.level1_short_pair_distance_ratio);
        configured_max_distance = ceil_to_cumulative_time_int_for_adaptive_bounds(scaled);
    }
    if (!config.level1_short_pair_acceptance || configured_max_distance == 0 ||
        result.residual_active_detectors.size() < 2) {
        return;
    }

    std::vector<uint64_t> residual = result.residual_active_detectors;
    std::sort(residual.begin(), residual.end());
    residual.erase(std::unique(residual.begin(), residual.end()), residual.end());
    const size_t n = residual.size();
    if (n < 2) {
        result.residual_active_detectors = std::move(residual);
        return;
    }

    const cumulative_time_int max_distance = configured_max_distance;
    const cumulative_time_int external_guard = config.level1_short_pair_external_guard_explicit
        ? config.level1_short_pair_external_guard
        : max_distance;
    const cumulative_time_int distance_cap = std::max(max_distance, external_guard);

    std::vector<std::tuple<cumulative_time_int, size_t, size_t>> candidates;
    candidates.reserve(n);
    std::vector<cumulative_time_int> nearest(n, GROWING_ONLY_INF_DISTANCE);
    std::vector<size_t> nearest_count(n, 0);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            cumulative_time_int distance = shortest_detector_distance_for_level1_short_pair(
                graph, graph_cache, residual[i], residual[j], distance_cap);
            if (distance >= GROWING_ONLY_INF_DISTANCE) {
                continue;
            }
            if (distance < nearest[i]) {
                nearest[i] = distance;
                nearest_count[i] = 1;
            } else if (distance == nearest[i]) {
                nearest_count[i]++;
            }
            if (distance < nearest[j]) {
                nearest[j] = distance;
                nearest_count[j] = 1;
            } else if (distance == nearest[j]) {
                nearest_count[j]++;
            }
            if (distance <= max_distance && distance <= bounds.diameter_bound) {
                candidates.push_back({distance, i, j});
            }
        }
    }
    if (candidates.empty()) {
        result.residual_active_detectors = std::move(residual);
        return;
    }
    std::sort(candidates.begin(), candidates.end());

    auto passes_external_guard = [&](size_t i, size_t j) {
        if (external_guard == 0) {
            return true;
        }
        for (uint64_t other : original_active_detectors) {
            if (other == residual[i] || other == residual[j]) {
                continue;
            }
            cumulative_time_int di = shortest_detector_distance_for_level1_short_pair(
                graph, graph_cache, residual[i], other, external_guard);
            if (di <= external_guard) {
                return false;
            }
            cumulative_time_int dj = shortest_detector_distance_for_level1_short_pair(
                graph, graph_cache, residual[j], other, external_guard);
            if (dj <= external_guard) {
                return false;
            }
        }
        return true;
    };

    std::vector<uint8_t> used(n, 0);
    std::vector<uint64_t> next_residual;
    size_t added_clusters = 0;
    size_t added_active = 0;
    for (const auto& candidate : candidates) {
        cumulative_time_int distance;
        size_t i;
        size_t j;
        std::tie(distance, i, j) = candidate;
        if (used[i] || used[j]) {
            continue;
        }
        if (config.level1_short_pair_require_mutual_nearest &&
            !(nearest[i] == distance && nearest[j] == distance)) {
            continue;
        }
        if (!passes_external_guard(i, j)) {
            continue;
        }

        used[i] = 1;
        used[j] = 1;
        GrowingOnlyAcceptedCluster cluster;
        cluster.active_detectors = {residual[i], residual[j]};
        std::sort(cluster.active_detectors.begin(), cluster.active_detectors.end());
        cluster.diameter = distance;
        cluster.max_internal_pair_distance = distance;
        cluster.max_internal_pair_collision_time = static_cast<cumulative_time_int>(distance / 2);
        cluster.ready_time = cluster.max_internal_pair_collision_time;
        cluster.first_external_collision_time = GROWING_ONLY_INF_DISTANCE;
        if (external_guard != 0) {
            cluster.first_external_collision_time = static_cast<cumulative_time_int>(external_guard / 2 + 1);
        }
        cluster.buffer = cluster.first_external_collision_time >= GROWING_ONLY_INF_DISTANCE
                             ? GROWING_ONLY_INF_DISTANCE
                             : (cluster.first_external_collision_time > cluster.ready_time
                                    ? cluster.first_external_collision_time - cluster.ready_time
                                    : 0);
        cluster.nearest_boundary_match_distance = GROWING_ONLY_INF_DISTANCE;
        cluster.parity_can_stop_locally = true;
        cluster.forced_by_max_level = false;
        cluster.internal_collisions_recorded = 1;
        cluster.external_collisions_recorded = 0;
        result.accepted_clusters.push_back(std::move(cluster));
        added_clusters++;
        added_active += 2;
    }

    for (size_t i = 0; i < n; i++) {
        if (!used[i]) {
            next_residual.push_back(residual[i]);
        }
    }
    if (added_clusters == 0) {
        result.residual_active_detectors = std::move(residual);
        return;
    }
    result.residual_active_detectors = std::move(next_residual);
    result.stats.removed_or_hidden_clusters += added_clusters;
    result.stats.max_component_active_size = std::max<size_t>(result.stats.max_component_active_size, 2);
    if (result.stats.max_residual_active_size < result.residual_active_detectors.size()) {
        result.stats.max_residual_active_size = result.residual_active_detectors.size();
    }
    // Reuse the deferral counters as generic adaptive-residual-change counters
    // would be misleading, so only expose this via trace output below.
    (void)added_active;
}

void apply_adaptive_heavy_cluster_deferral(
    GrowingOnlyClustererResult& result,
    const ProcessingClusterConfig& config) {
    if (!config.adaptive_defer_heavy_clusters || result.accepted_clusters.empty()) {
        return;
    }
    std::vector<GrowingOnlyAcceptedCluster> kept;
    kept.reserve(result.accepted_clusters.size());
    size_t deferred_clusters = 0;
    size_t deferred_active = 0;
    for (auto& cluster : result.accepted_clusters) {
        if (adaptive_cluster_should_be_deferred(cluster, config)) {
            deferred_clusters++;
            deferred_active += cluster.active_detectors.size();
            for (auto detector : cluster.active_detectors) {
                result.residual_active_detectors.push_back(detector);
            }
        } else {
            kept.push_back(cluster);
        }
    }
    if (deferred_clusters == 0) {
        return;
    }
    result.accepted_clusters = std::move(kept);
    std::sort(result.residual_active_detectors.begin(), result.residual_active_detectors.end());
    result.residual_active_detectors.erase(
        std::unique(result.residual_active_detectors.begin(), result.residual_active_detectors.end()),
        result.residual_active_detectors.end());
    result.stats.adaptive_deferred_clusters += deferred_clusters;
    result.stats.adaptive_deferred_active_detectors += deferred_active;
}

GrowingOnlyClustererResult run_growing_only_clustering_for_adaptive_trial(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& active_detectors,
    size_t level,
    const ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache) {
    auto result = run_growing_only_clustering(graph, active_detectors, level, config, graph_cache);
    apply_adaptive_heavy_cluster_deferral(result, config);
    return result;
}

struct AdaptiveMaxActiveRolloutCost {
    size_t worst_active = 0;
    uint64_t worst_events = 0;
    uint64_t ideal_events = 0;
    std::vector<uint64_t> level_critical_events;
    size_t nonempty_event_levels = 0;
    size_t residual_after = 0;
    size_t accepted_active_total = 0;
    size_t accepted_clusters_total = 0;
    cumulative_time_int final_stop_time = 0;
    bool have_completed_level = false;
};

double adaptive_weighted_ideal_events_from_levels(
    const std::vector<uint64_t>& level_events,
    double weight_power) {
    if (level_events.empty()) {
        return 0.0;
    }
    uint64_t max_event = 0;
    for (uint64_t e : level_events) {
        max_event = std::max(max_event, e);
    }
    if (max_event == 0) {
        return 0.0;
    }
    if (!(weight_power > 0.0)) {
        double total = 0.0;
        for (uint64_t e : level_events) {
            total += static_cast<double>(e);
        }
        return total;
    }
    const double denom = static_cast<double>(max_event);
    double total = 0.0;
    for (uint64_t e : level_events) {
        if (e == 0) {
            continue;
        }
        const double ratio = static_cast<double>(e) / denom;
        total += static_cast<double>(e) * std::pow(ratio, weight_power);
    }
    return total;
}

GrowingOnlyClustererResult choose_and_run_breakpoint_adaptive_growing_only_clustering(
    const MatchingGraph& graph,
    std::vector<uint64_t> active_detectors,
    size_t level,
    ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache,
    size_t lookahead_depth_override);

struct AdaptiveExactLevelAccounting {
    uint64_t prior_critical_events = 0;
    cumulative_time_int prior_stop_time = 0;
    bool have_prior_completed_level = false;
    uint64_t next_critical_events = 0;
    cumulative_time_int next_stop_time = 0;
    bool have_completed_level = false;
    uint64_t level_max_additional_events = 0;
    uint64_t level_max_raw_events = 0;
    size_t max_active = 0;
    size_t cluster_count = 0;
    size_t active_total = 0;
};

AdaptiveExactLevelAccounting adaptive_exact_account_detector_sets_for_scheduler_metric(
    const std::vector<std::vector<uint64_t>>& detector_sets,
    const ProcessingClusterConfig& config) {
    AdaptiveExactLevelAccounting out;
    out.prior_critical_events = config.adaptive_exact_ideal_prior_critical_events;
    out.prior_stop_time = config.adaptive_exact_ideal_prior_stop_time;
    out.have_prior_completed_level = config.adaptive_exact_ideal_have_prior_completed_level;
    out.next_critical_events = out.prior_critical_events;
    out.next_stop_time = out.prior_stop_time;
    out.have_completed_level = out.have_prior_completed_level;
    if (detector_sets.empty()) {
        return out;
    }
    const cumulative_time_int cutoff_time = out.have_prior_completed_level
        ? out.prior_stop_time
        : std::numeric_limits<cumulative_time_int>::min();
    cumulative_time_int level_max_stop_time = out.prior_stop_time;
    bool any_nonempty = false;
    for (const auto& detectors : detector_sets) {
        if (detectors.empty()) {
            continue;
        }
        any_nonempty = true;
        out.cluster_count++;
        out.active_total += detectors.size();
        out.max_active = std::max(out.max_active, detectors.size());
        const AdaptiveEventTiming timing = adaptive_event_timing_for_detectors(config, detectors, cutoff_time);
        out.level_max_additional_events = std::max(out.level_max_additional_events, timing.events_after_cutoff);
        out.level_max_raw_events = std::max(out.level_max_raw_events, timing.raw_event_count);
        level_max_stop_time = std::max(level_max_stop_time, timing.stop_algorithmic_time);
    }
    if (any_nonempty) {
        out.next_critical_events = out.prior_critical_events + out.level_max_additional_events;
        out.next_stop_time = std::max(out.prior_stop_time, level_max_stop_time);
        out.have_completed_level = true;
    }
    return out;
}

std::vector<std::vector<uint64_t>> adaptive_detector_sets_from_result(
    const GrowingOnlyClustererResult& result) {
    std::vector<std::vector<uint64_t>> detector_sets;
    detector_sets.reserve(result.accepted_clusters.size());
    for (const auto& cluster : result.accepted_clusters) {
        detector_sets.push_back(cluster.active_detectors);
    }
    return detector_sets;
}

ProcessingClusterConfig adaptive_config_with_exact_prior(
    ProcessingClusterConfig config,
    uint64_t critical_events,
    cumulative_time_int stop_time,
    bool have_completed_level) {
    config.adaptive_exact_ideal_prior_critical_events = critical_events;
    config.adaptive_exact_ideal_prior_stop_time = stop_time;
    config.adaptive_exact_ideal_have_prior_completed_level = have_completed_level;
    return config;
}

AdaptiveMaxActiveRolloutCost adaptive_lookahead_max_active_rollout_cost(
    const MatchingGraph& graph,
    std::vector<uint64_t> residual,
    size_t level,
    ProcessingClusterConfig config,
    const ProcessingClusterGraphCache* graph_cache,
    size_t remaining_depth) {
    std::sort(residual.begin(), residual.end());
    residual.erase(std::unique(residual.begin(), residual.end()), residual.end());
    const bool use_exact_scheduler_metric =
        config.adaptive_objective_minimize_ideal_events &&
        static_cast<bool>(config.adaptive_event_timing_callback);
    if (residual.empty()) {
        AdaptiveMaxActiveRolloutCost cost;
        if (use_exact_scheduler_metric) {
            cost.ideal_events = config.adaptive_exact_ideal_prior_critical_events;
            cost.final_stop_time = config.adaptive_exact_ideal_prior_stop_time;
            cost.have_completed_level = config.adaptive_exact_ideal_have_prior_completed_level;
        }
        return cost;
    }
    const size_t max_level = processing_cluster_effective_max_level(config);

    if (use_exact_scheduler_metric) {
        if (level > max_level) {
            std::vector<std::vector<uint64_t>> root_sets;
            root_sets.push_back(residual);
            const auto accounted = adaptive_exact_account_detector_sets_for_scheduler_metric(root_sets, config);
            AdaptiveMaxActiveRolloutCost cost;
            cost.worst_active = residual.size();
            cost.worst_events = accounted.level_max_raw_events;
            cost.ideal_events = accounted.next_critical_events;
            if (accounted.level_max_additional_events > 0) {
                cost.level_critical_events.push_back(accounted.level_max_additional_events);
            }
            cost.nonempty_event_levels = accounted.cluster_count > 0 ? 1 : 0;
            cost.residual_after = residual.size();
            cost.accepted_active_total = accounted.active_total;
            cost.accepted_clusters_total = accounted.cluster_count;
            cost.final_stop_time = accounted.next_stop_time;
            cost.have_completed_level = accounted.have_completed_level;
            return cost;
        }

        config.adaptive_trace_bounds = false;
        const size_t next_depth = remaining_depth == 0 ? 0 : remaining_depth - 1;
        auto result = choose_and_run_breakpoint_adaptive_growing_only_clustering(
            graph,
            std::move(residual),
            level,
            config,
            graph_cache,
            next_depth);

        const auto detector_sets = adaptive_detector_sets_from_result(result);
        const auto accounted = adaptive_exact_account_detector_sets_for_scheduler_metric(detector_sets, config);
        ProcessingClusterConfig tail_config = adaptive_config_with_exact_prior(
            config,
            accounted.next_critical_events,
            accounted.next_stop_time,
            accounted.have_completed_level);
        AdaptiveMaxActiveRolloutCost tail = adaptive_lookahead_max_active_rollout_cost(
            graph,
            result.residual_active_detectors,
            level + 1,
            tail_config,
            graph_cache,
            next_depth);

        AdaptiveMaxActiveRolloutCost cost;
        cost.worst_active = std::max(accounted.max_active, tail.worst_active);
        cost.worst_events = std::max(accounted.level_max_raw_events, tail.worst_events);
        cost.ideal_events = tail.ideal_events;
        if (accounted.level_max_additional_events > 0) {
            cost.level_critical_events.push_back(accounted.level_max_additional_events);
        }
        cost.level_critical_events.insert(
            cost.level_critical_events.end(),
            tail.level_critical_events.begin(),
            tail.level_critical_events.end());
        cost.nonempty_event_levels = tail.nonempty_event_levels + (accounted.cluster_count > 0 ? 1 : 0);
        cost.residual_after = tail.residual_after;
        cost.accepted_active_total = accounted.active_total + tail.accepted_active_total;
        cost.accepted_clusters_total = accounted.cluster_count + tail.accepted_clusters_total;
        cost.final_stop_time = tail.final_stop_time;
        cost.have_completed_level = tail.have_completed_level;
        return cost;
    }

    if (level > max_level) {
        const uint64_t residual_events = adaptive_event_count_proxy_for_detectors(config, residual);
        AdaptiveMaxActiveRolloutCost cost;
        cost.worst_active = residual.size();
        cost.worst_events = residual_events;
        cost.ideal_events = residual_events;
        if (residual_events > 0) {
            cost.level_critical_events.push_back(residual_events);
        }
        cost.nonempty_event_levels = residual_events > 0 ? 1 : 0;
        cost.residual_after = residual.size();
        return cost;
    }
    if (remaining_depth == 0) {
        const uint64_t residual_events = adaptive_event_count_proxy_for_detectors(config, residual);
        AdaptiveMaxActiveRolloutCost cost;
        cost.worst_active = residual.size();
        cost.worst_events = residual_events;
        cost.ideal_events = residual_events;
        if (residual_events > 0) {
            cost.level_critical_events.push_back(residual_events);
        }
        cost.nonempty_event_levels = residual_events > 0 ? 1 : 0;
        cost.residual_after = residual.size();
        return cost;
    }

    config.adaptive_trace_bounds = false;
    auto result = choose_and_run_breakpoint_adaptive_growing_only_clustering(
        graph,
        std::move(residual),
        level,
        config,
        graph_cache,
        remaining_depth - 1);

    AdaptiveMaxActiveRolloutCost tail = adaptive_lookahead_max_active_rollout_cost(
        graph,
        result.residual_active_detectors,
        level + 1,
        config,
        graph_cache,
        remaining_depth - 1);

    AdaptiveMaxActiveRolloutCost cost;
    const uint64_t current_max_events = max_accepted_cluster_events_for_adaptive_result(result, config);
    cost.worst_active = std::max(max_accepted_cluster_active_for_adaptive_result(result), tail.worst_active);
    cost.worst_events = std::max(current_max_events, tail.worst_events);
    cost.ideal_events = std::max(current_max_events, tail.ideal_events);
    if (current_max_events > 0) {
        cost.level_critical_events.push_back(current_max_events);
    }
    cost.level_critical_events.insert(
        cost.level_critical_events.end(),
        tail.level_critical_events.begin(),
        tail.level_critical_events.end());
    cost.nonempty_event_levels = tail.nonempty_event_levels + (current_max_events > 0 ? 1 : 0);
    cost.residual_after = tail.residual_after;
    cost.accepted_active_total = accepted_active_count_for_adaptive_result(result) + tail.accepted_active_total;
    cost.accepted_clusters_total = result.accepted_clusters.size() + tail.accepted_clusters_total;
    return cost;
}

struct AdaptiveMaxActiveCandidateScore {
    size_t worst_active = std::numeric_limits<size_t>::max();
    uint64_t worst_events = std::numeric_limits<uint64_t>::max();
    uint64_t ideal_events = std::numeric_limits<uint64_t>::max();
    double weighted_ideal_events = std::numeric_limits<double>::infinity();
    size_t nonempty_event_levels = std::numeric_limits<size_t>::max();
    size_t residual_after_lookahead = std::numeric_limits<size_t>::max();
    size_t accepted_active_total = 0;
    size_t current_max_active = std::numeric_limits<size_t>::max();
    uint64_t current_max_events = std::numeric_limits<uint64_t>::max();
    size_t accepted_clusters_total = 0;
    size_t current_accepted_clusters = 0;
    cumulative_time_int diameter_bound = GROWING_ONLY_INF_DISTANCE;
    cumulative_time_int buffer_bound = GROWING_ONLY_INF_DISTANCE;
};

AdaptiveMaxActiveCandidateScore adaptive_max_active_score_for_candidate(
    const MatchingGraph& graph,
    const GrowingOnlyClustererResult& current,
    const ProcessingClusterBounds& current_bounds,
    size_t level,
    const ProcessingClusterConfig& config_after_current_bounds,
    const ProcessingClusterGraphCache* graph_cache,
    size_t lookahead_depth) {
    const bool use_exact_scheduler_metric =
        config_after_current_bounds.adaptive_objective_minimize_ideal_events &&
        static_cast<bool>(config_after_current_bounds.adaptive_event_timing_callback);

    AdaptiveMaxActiveCandidateScore score;
    score.current_max_active = max_accepted_cluster_active_for_adaptive_result(current);
    score.current_accepted_clusters = current.accepted_clusters.size();
    score.diameter_bound = current_bounds.diameter_bound;
    score.buffer_bound = current_bounds.buffer_bound;

    if (use_exact_scheduler_metric) {
        const auto detector_sets = adaptive_detector_sets_from_result(current);
        const auto accounted = adaptive_exact_account_detector_sets_for_scheduler_metric(
            detector_sets,
            config_after_current_bounds);
        ProcessingClusterConfig future_config = adaptive_config_with_exact_prior(
            config_after_current_bounds,
            accounted.next_critical_events,
            accounted.next_stop_time,
            accounted.have_completed_level);
        AdaptiveMaxActiveRolloutCost future = adaptive_lookahead_max_active_rollout_cost(
            graph,
            current.residual_active_detectors,
            level + 1,
            future_config,
            graph_cache,
            lookahead_depth);

        score.current_max_events = accounted.level_max_raw_events;
        score.worst_active = std::max(accounted.max_active, future.worst_active);
        score.worst_events = std::max(accounted.level_max_raw_events, future.worst_events);
        score.ideal_events = future.ideal_events;
        score.weighted_ideal_events = static_cast<double>(score.ideal_events);
        score.nonempty_event_levels = future.nonempty_event_levels + (accounted.cluster_count > 0 ? 1 : 0);
        score.residual_after_lookahead = future.residual_after;
        score.accepted_active_total = accounted.active_total + future.accepted_active_total;
        score.accepted_clusters_total = accounted.cluster_count + future.accepted_clusters_total;
        return score;
    }

    AdaptiveMaxActiveRolloutCost future = adaptive_lookahead_max_active_rollout_cost(
        graph,
        current.residual_active_detectors,
        level + 1,
        config_after_current_bounds,
        graph_cache,
        lookahead_depth);

    score.current_max_events = max_accepted_cluster_events_for_adaptive_result(current, config_after_current_bounds);
    score.worst_active = std::max(score.current_max_active, future.worst_active);
    score.worst_events = std::max(score.current_max_events, future.worst_events);
    score.ideal_events = std::max(score.current_max_events, future.ideal_events);
    std::vector<uint64_t> candidate_level_events;
    if (score.current_max_events > 0) {
        candidate_level_events.push_back(score.current_max_events);
    }
    candidate_level_events.insert(
        candidate_level_events.end(),
        future.level_critical_events.begin(),
        future.level_critical_events.end());
    (void)candidate_level_events;
    score.weighted_ideal_events = static_cast<double>(score.ideal_events);
    score.nonempty_event_levels = future.nonempty_event_levels + (score.current_max_events > 0 ? 1 : 0);
    score.residual_after_lookahead = future.residual_after;
    score.accepted_active_total = accepted_active_count_for_adaptive_result(current) + future.accepted_active_total;
    score.accepted_clusters_total = current.accepted_clusters.size() + future.accepted_clusters_total;
    return score;
}

bool adaptive_max_active_score_is_better(
    const AdaptiveMaxActiveCandidateScore& candidate,
    const AdaptiveMaxActiveCandidateScore& best,
    bool have_best,
    bool prefer_events,
    bool prefer_ideal_events) {
    if (!have_best) {
        return true;
    }
    if (prefer_ideal_events) {
        if (candidate.ideal_events != best.ideal_events) {
            return candidate.ideal_events < best.ideal_events;
        }
        if (candidate.worst_events != best.worst_events) {
            return candidate.worst_events < best.worst_events;
        }
        if (candidate.residual_after_lookahead != best.residual_after_lookahead) {
            return candidate.residual_after_lookahead < best.residual_after_lookahead;
        }
        if (candidate.accepted_active_total != best.accepted_active_total) {
            return candidate.accepted_active_total > best.accepted_active_total;
        }
        if (candidate.current_max_events != best.current_max_events) {
            return candidate.current_max_events < best.current_max_events;
        }
    } else if (prefer_events) {
        if (candidate.worst_events != best.worst_events) {
            return candidate.worst_events < best.worst_events;
        }
        if (candidate.residual_after_lookahead != best.residual_after_lookahead) {
            return candidate.residual_after_lookahead < best.residual_after_lookahead;
        }
        if (candidate.accepted_active_total != best.accepted_active_total) {
            return candidate.accepted_active_total > best.accepted_active_total;
        }
        if (candidate.current_max_events != best.current_max_events) {
            return candidate.current_max_events < best.current_max_events;
        }
    }
    if (candidate.worst_active != best.worst_active) {
        return candidate.worst_active < best.worst_active;
    }
    if (candidate.residual_after_lookahead != best.residual_after_lookahead) {
        return candidate.residual_after_lookahead < best.residual_after_lookahead;
    }
    if (candidate.accepted_active_total != best.accepted_active_total) {
        return candidate.accepted_active_total > best.accepted_active_total;
    }
    if (candidate.current_max_active != best.current_max_active) {
        return candidate.current_max_active < best.current_max_active;
    }
    if (candidate.accepted_clusters_total != best.accepted_clusters_total) {
        return candidate.accepted_clusters_total > best.accepted_clusters_total;
    }
    if (candidate.current_accepted_clusters != best.current_accepted_clusters) {
        return candidate.current_accepted_clusters > best.current_accepted_clusters;
    }
    if (candidate.diameter_bound != best.diameter_bound) {
        return candidate.diameter_bound < best.diameter_bound;
    }
    return candidate.buffer_bound < best.buffer_bound;
}

bool adaptive_result_score_is_better(
    const GrowingOnlyClustererResult& candidate,
    const ProcessingClusterBounds& candidate_bounds,
    const GrowingOnlyClustererResult& best,
    const ProcessingClusterBounds& best_bounds,
    bool have_best) {
    if (!have_best) {
        return true;
    }
    const size_t c_clusters = candidate.accepted_clusters.size();
    const size_t b_clusters = best.accepted_clusters.size();
    if (c_clusters != b_clusters) {
        return c_clusters > b_clusters;
    }
    const size_t c_active = accepted_active_count_for_adaptive_result(candidate);
    const size_t b_active = accepted_active_count_for_adaptive_result(best);
    if (c_active != b_active) {
        return c_active > b_active;
    }
    if (candidate.residual_active_detectors.size() != best.residual_active_detectors.size()) {
        return candidate.residual_active_detectors.size() < best.residual_active_detectors.size();
    }
    if (candidate_bounds.diameter_bound != best_bounds.diameter_bound) {
        return candidate_bounds.diameter_bound < best_bounds.diameter_bound;
    }
    return candidate_bounds.buffer_bound < best_bounds.buffer_bound;
}

cumulative_time_int adaptive_required_diameter_for_accepted_cluster(
    const GrowingOnlyAcceptedCluster& cluster) {
    cumulative_time_int required = cluster.diameter;
    if ((cluster.active_detectors.size() & 1) != 0) {
        required = std::max(required, cluster.nearest_boundary_match_distance);
    }
    return required;
}

cumulative_time_int adaptive_required_diameter_for_result(
    const GrowingOnlyClustererResult& result,
    cumulative_time_int lower_bound) {
    cumulative_time_int required = lower_bound;
    for (const auto& cluster : result.accepted_clusters) {
        required = std::max(required, adaptive_required_diameter_for_accepted_cluster(cluster));
    }
    return required;
}


std::vector<cumulative_time_int> adaptive_level1_diameter_candidates_from_active_detectors(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& active_detectors,
    const ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache,
    const ProcessingClusterBounds& fixed_bounds) {
    std::vector<cumulative_time_int> values;
    auto add_value = [&](cumulative_time_int value) {
        add_positive_candidate_distance(values, value);
    };

    const cumulative_time_int fixed_d = fixed_bounds.diameter_bound;
    add_value(fixed_d);

    cumulative_time_int max_d = config.adaptive_level1_max_diameter;
    if (max_d == 0) {
        long double ratio = config.adaptive_level1_max_diameter_ratio;
        if (!(ratio > 0.0) || !std::isfinite(static_cast<double>(ratio))) {
            ratio = 3.0L;
        }
        long double scaled = static_cast<long double>(fixed_d) * ratio;
        max_d = ceil_to_cumulative_time_int_for_adaptive_bounds(scaled);
    }
    if (max_d == 0 || max_d >= GROWING_ONLY_INF_DISTANCE) {
        max_d = fixed_d;
    }
    max_d = std::max(max_d, fixed_d);

    if (graph_cache != nullptr && !graph_cache->nearest_boundary_match_distance_by_vertex.empty()) {
        for (auto detector : active_detectors) {
            if (detector < graph_cache->nearest_boundary_match_distance_by_vertex.size()) {
                cumulative_time_int boundary = graph_cache->nearest_boundary_match_distance_by_vertex[detector];
                if (boundary <= max_d) {
                    add_value(boundary);
                }
            }
        }
    } else {
        (void)graph;
    }

    const size_t n = active_detectors.size();
    if (graph_cache != nullptr && graph_cache->has_all_pairs_distances) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = i + 1; j < n; j++) {
                cumulative_time_int distance = cached_all_pairs_distance_for_adaptive_bounds(
                    graph_cache, active_detectors[i], active_detectors[j]);
                if (distance <= max_d) {
                    add_value(distance);
                }
            }
        }
    } else {
        // No all-pairs cache: run one bounded Dijkstra per active detector,
        // instead of one Dijkstra per active pair.  This keeps level-1
        // re-optimization usable when all-pairs preprocessing is disabled.
        std::vector<size_t> active_indices;
        active_indices.reserve(n);
        for (auto detector : active_detectors) {
            active_indices.push_back(static_cast<size_t>(detector));
        }
        for (size_t i = 0; i < n; i++) {
            const size_t source = active_indices[i];
            if (source >= graph.nodes.size()) {
                continue;
            }
            std::vector<cumulative_time_int> distances(graph.nodes.size(), GROWING_ONLY_INF_DISTANCE);
            using QueueEntry = std::pair<cumulative_time_int, size_t>;
            std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
            distances[source] = 0;
            queue.push({0, source});
            while (!queue.empty()) {
                auto [distance, node_index] = queue.top();
                queue.pop();
                if (distance != distances[node_index]) {
                    continue;
                }
                if (distance > max_d) {
                    break;
                }
                const auto& node = graph.nodes[node_index];
                for (size_t edge_index = 0; edge_index < node.neighbors.size(); edge_index++) {
                    const auto* neighbor = node.neighbors[edge_index];
                    if (neighbor == nullptr) {
                        continue;
                    }
                    size_t next_index = static_cast<size_t>(neighbor - graph.nodes.data());
                    if (next_index >= graph.nodes.size()) {
                        continue;
                    }
                    cumulative_time_int weight = static_cast<cumulative_time_int>(node.neighbor_weights[edge_index]);
                    if (distance >= GROWING_ONLY_INF_DISTANCE - weight) {
                        continue;
                    }
                    cumulative_time_int next_distance = distance + weight;
                    if (next_distance < distances[next_index] && next_distance <= max_d) {
                        distances[next_index] = next_distance;
                        queue.push({next_distance, next_index});
                    }
                }
            }
            for (size_t j = i + 1; j < n; j++) {
                size_t destination = active_indices[j];
                if (destination < distances.size() && distances[destination] <= max_d) {
                    add_value(distances[destination]);
                }
            }
        }
    }

    size_t limit = config.adaptive_level1_candidate_limit;
    if (limit == 0) {
        limit = std::max<size_t>(1, config.adaptive_candidate_limit);
    }
    return unique_limited_adaptive_values(std::move(values), limit);
}

std::vector<cumulative_time_int> adaptive_breakpoint_buffer_candidates_from_residual(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& active_detectors,
    const ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache,
    size_t level) {
    std::vector<cumulative_time_int> values;
    auto add_value = [&](cumulative_time_int value) {
        add_positive_candidate_distance(values, value);
    };

    if (level > 0 && level <= config.explicit_bounds_by_level.size()) {
        add_value(config.explicit_bounds_by_level[level - 1].buffer_bound);
    }
    if (level > 1 && level - 2 < config.explicit_bounds_by_level.size()) {
        const auto& previous = config.explicit_bounds_by_level[level - 2];
        add_value(previous.buffer_bound);
        add_value(previous.buffer_bound + 1);
        add_value(adaptive_required_buffer_for_diameter(
            config, level, adaptive_legacy_gap_lower_bound(config, level)));
    }

    // Singletons and odd components are controlled by their boundary-inclusive
    // required diameter.  Add the exact certificate buffer required for those
    // distances so a sparse shot with no close detector-detector pairs still has
    // a meaningful residual-derived breakpoint candidate.
    if (graph_cache != nullptr && !graph_cache->nearest_boundary_match_distance_by_vertex.empty()) {
        for (auto detector : active_detectors) {
            if (detector < graph_cache->nearest_boundary_match_distance_by_vertex.size()) {
                cumulative_time_int boundary = graph_cache->nearest_boundary_match_distance_by_vertex[detector];
                add_value(boundary);
                add_value(adaptive_required_buffer_for_diameter(
                    config,
                    level,
                    std::max(boundary, adaptive_legacy_gap_lower_bound(config, level))));
            }
        }
    } else {
        (void)graph;
    }

    const size_t active_count = active_detectors.size();
    for (size_t i = 0; i < active_count; i++) {
        for (size_t j = i + 1; j < active_count; j++) {
            cumulative_time_int distance = cached_all_pairs_distance_for_adaptive_bounds(
                graph_cache, active_detectors[i], active_detectors[j]);
            if (distance <= 0 || distance >= GROWING_ONLY_INF_DISTANCE) {
                continue;
            }
            // b=distance is immediately after this merge; b=distance-1 is the
            // last integer value before it.  These two values enumerate the
            // threshold-component breakpoints without scanning the whole range.
            add_value(distance);
            if (distance > 1) {
                add_value(distance - 1);
            }
            add_value(adaptive_required_buffer_for_diameter(
                config,
                level,
                std::max(distance, adaptive_legacy_gap_lower_bound(config, level))));
        }
    }

    if (values.empty()) {
        add_value(adaptive_required_buffer_for_diameter(
            config, level, adaptive_legacy_gap_lower_bound(config, level)));
        if (level > 0 && level <= config.explicit_bounds_by_level.size()) {
            add_value(config.explicit_bounds_by_level[level - 1].buffer_bound);
        }
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

GrowingOnlyClustererResult choose_and_run_breakpoint_adaptive_growing_only_clustering(
    const MatchingGraph& graph,
    std::vector<uint64_t> active_detectors,
    size_t level,
    ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache,
    size_t lookahead_depth_override) {
    ensure_adaptive_explicit_bounds_vector(config);
    std::sort(active_detectors.begin(), active_detectors.end());
    active_detectors.erase(std::unique(active_detectors.begin(), active_detectors.end()), active_detectors.end());

    // Level 1 can either stay fixed (legacy behavior) or be re-optimized from
    // scratch.  In re-optimization mode we do not post-hoc accept short pairs;
    // instead each candidate d_1 recomputes a certified b_1 and reruns ordinary
    // level-1 clustering with that actual (d_1,b_1).  This keeps additional
    // o-x-o / o-x-x-o absorption inside the normal certificate.
    if (level <= 1) {
        ProcessingClusterBounds fixed_bounds = config.explicit_bounds_by_level[level - 1];
        adaptive_fill_certificate_metadata(config, level, fixed_bounds);

        if (config.adaptive_level1_reoptimize_bounds) {
            const auto d_values = adaptive_level1_diameter_candidates_from_active_detectors(
                graph, active_detectors, config, graph_cache, fixed_bounds);
            GrowingOnlyClustererResult best_result;
            ProcessingClusterBounds best_bounds{};
            AdaptiveMaxActiveCandidateScore best_max_active_score;
            bool have_best = false;
            size_t trials = 0;

            for (auto d : d_values) {
                if (!(d > 0) || d >= GROWING_ONLY_INF_DISTANCE) {
                    continue;
                }
                ProcessingClusterBounds candidate_bounds{d, adaptive_required_buffer_for_diameter(config, level, d)};
                adaptive_fill_certificate_metadata(config, level, candidate_bounds);
                if (!adaptive_bounds_satisfy_optional_constraints(config, level, candidate_bounds)) {
                    continue;
                }
                auto candidate_config = adaptive_trial_config_with_bounds(config, level, candidate_bounds);
                auto candidate_result = run_growing_only_clustering_for_adaptive_trial(
                    graph, active_detectors, level, candidate_config, graph_cache);
                trials++;

                bool take_candidate = false;
                if (config.adaptive_objective_minimize_max_active || config.adaptive_objective_minimize_max_events || config.adaptive_objective_minimize_ideal_events) {
                    ProcessingClusterConfig score_config = config;
                    ensure_adaptive_explicit_bounds_vector(score_config);
                    ProcessingClusterBounds score_bounds = candidate_bounds;
                    adaptive_fill_certificate_metadata(score_config, level, score_bounds);
                    score_config.explicit_bounds_by_level[level - 1] = score_bounds;
                    auto score = adaptive_max_active_score_for_candidate(
                        graph,
                        candidate_result,
                        score_bounds,
                        level,
                        score_config,
                        graph_cache,
                        lookahead_depth_override);
                    take_candidate = adaptive_max_active_score_is_better(
                        score,
                        best_max_active_score,
                        have_best,
                        config.adaptive_objective_minimize_max_events,
                        config.adaptive_objective_minimize_ideal_events);
                    if (take_candidate) {
                        best_max_active_score = score;
                    }
                } else {
                    take_candidate = adaptive_result_score_is_better(
                        candidate_result, candidate_bounds, best_result, best_bounds, have_best);
                }
                if (take_candidate) {
                    best_result = std::move(candidate_result);
                    best_bounds = candidate_bounds;
                    have_best = true;
                }
            }

            if (!have_best) {
                best_bounds = fixed_bounds;
                adaptive_raise_buffer_to_optional_lower_bounds(config, level, best_bounds);
                adaptive_fill_certificate_metadata(config, level, best_bounds);
                auto fallback_config = adaptive_trial_config_with_bounds(config, level, best_bounds);
                best_result = run_growing_only_clustering_for_adaptive_trial(
                    graph, active_detectors, level, fallback_config, graph_cache);
                if (config.adaptive_objective_minimize_max_active || config.adaptive_objective_minimize_max_events || config.adaptive_objective_minimize_ideal_events) {
                    ProcessingClusterConfig score_config = config;
                    ensure_adaptive_explicit_bounds_vector(score_config);
                    ProcessingClusterBounds score_bounds = best_bounds;
                    adaptive_fill_certificate_metadata(score_config, level, score_bounds);
                    score_config.explicit_bounds_by_level[level - 1] = score_bounds;
                    best_max_active_score = adaptive_max_active_score_for_candidate(
                        graph,
                        best_result,
                        score_bounds,
                        level,
                        score_config,
                        graph_cache,
                        lookahead_depth_override);
                }
            }

            config.explicit_bounds_by_level[level - 1] = best_bounds;
            best_result.bounds = best_bounds;
            if (config.adaptive_trace_bounds) {
                std::cerr << "ADAPTIVE_CLUSTER_BOUNDS_SELECTED"
                          << " strategy=breakpoint"
                          << " level=" << level
                          << " input_active=" << active_detectors.size()
                          << " d=" << best_bounds.diameter_bound
                          << " b=" << best_bounds.buffer_bound
                          << " accepted_clusters=" << best_result.accepted_clusters.size()
                          << " accepted_active=" << accepted_active_count_for_adaptive_result(best_result)
                          << " residual_active=" << best_result.residual_active_detectors.size()
                          << " trials=" << trials
                          << " fixed_level1=0"
                          << " level1_reopt=1"
                          << " level1_reopt_candidates=" << d_values.size()
                          << " level1_reopt_max_ratio=" << config.adaptive_level1_max_diameter_ratio
                          << " level1_reopt_max_d=" << config.adaptive_level1_max_diameter
                          << " objective=" << adaptive_objective_label(config)
                          << " lookahead=" << lookahead_depth_override
                          << " max_active=" << max_accepted_cluster_active_for_adaptive_result(best_result)
                          << " max_events=" << max_accepted_cluster_events_for_adaptive_result(best_result, config)
                          << " score_worst_active=" << ((config.adaptive_objective_minimize_max_active || config.adaptive_objective_minimize_max_events || config.adaptive_objective_minimize_ideal_events) ? best_max_active_score.worst_active : 0)
                          << " score_worst_events=" << ((config.adaptive_objective_minimize_max_events || config.adaptive_objective_minimize_ideal_events) ? best_max_active_score.worst_events : 0)
                          << " score_ideal_events=" << (config.adaptive_objective_minimize_ideal_events ? best_max_active_score.ideal_events : 0)
                          << " score_weighted_ideal_events=" << (config.adaptive_objective_minimize_ideal_events ? best_max_active_score.weighted_ideal_events : 0.0)
                          << " score_nonempty_event_levels=" << (config.adaptive_objective_minimize_ideal_events ? best_max_active_score.nonempty_event_levels : 0)
                          << " score_residual_after_lookahead=" << ((config.adaptive_objective_minimize_max_active || config.adaptive_objective_minimize_max_events || config.adaptive_objective_minimize_ideal_events) ? best_max_active_score.residual_after_lookahead : 0)
                          << " phi_cert_ok=" << (best_bounds.has_phi_buffer_certificate ? 1 : 0)
                          << " phi_required_b=" << best_bounds.phi_required_buffer_bound
                          << " ratio_required_b=" << best_bounds.ratio_required_buffer_bound
                          << " min_required_b=" << best_bounds.min_required_buffer_bound
                          << "\n";
            }
            return best_result;
        }

        config.explicit_bounds_by_level[level - 1] = fixed_bounds;
        auto fixed_config = adaptive_trial_config_with_bounds(config, level, fixed_bounds);
        auto fixed_result = run_growing_only_clustering_for_adaptive_trial(graph, active_detectors, level, fixed_config, graph_cache);
        const size_t accepted_before_short_pair = fixed_result.accepted_clusters.size();
        const size_t residual_before_short_pair = fixed_result.residual_active_detectors.size();
        apply_level1_short_pair_acceptance(
            graph, active_detectors, fixed_result, fixed_bounds, config, graph_cache);
        fixed_result.bounds = fixed_bounds;
        if (config.adaptive_trace_bounds) {
            std::cerr << "ADAPTIVE_CLUSTER_BOUNDS_SELECTED"
                      << " strategy=breakpoint"
                      << " level=" << level
                      << " input_active=" << active_detectors.size()
                      << " d=" << fixed_bounds.diameter_bound
                      << " b=" << fixed_bounds.buffer_bound
                      << " accepted_clusters=" << fixed_result.accepted_clusters.size()
                      << " accepted_active=" << accepted_active_count_for_adaptive_result(fixed_result)
                      << " residual_active=" << fixed_result.residual_active_detectors.size()
                      << " trials=0 fixed_level1=1"
                      << " level1_reopt=0"
                      << " level1_short_pair=" << (config.level1_short_pair_acceptance ? 1 : 0)
                      << " level1_short_pair_max_distance=" << (config.level1_short_pair_max_distance != 0 ? config.level1_short_pair_max_distance : static_cast<cumulative_time_int>(std::ceil(static_cast<long double>(fixed_bounds.diameter_bound) * static_cast<long double>(config.level1_short_pair_distance_ratio))))
                      << " level1_short_pair_distance_ratio=" << config.level1_short_pair_distance_ratio
                      << " level1_short_pair_external_guard=" << (config.level1_short_pair_external_guard_explicit ? config.level1_short_pair_external_guard : (config.level1_short_pair_max_distance != 0 ? config.level1_short_pair_max_distance : static_cast<cumulative_time_int>(std::ceil(static_cast<long double>(fixed_bounds.diameter_bound) * static_cast<long double>(config.level1_short_pair_distance_ratio)))))
                      << " level1_short_pair_added_clusters=" << (fixed_result.accepted_clusters.size() - accepted_before_short_pair)
                      << " level1_short_pair_residual_before=" << residual_before_short_pair
                      << "\n";
        }
        return fixed_result;
    }

    const cumulative_time_int d_lower = adaptive_legacy_gap_lower_bound(config, level);
    const auto b_values = adaptive_breakpoint_buffer_candidates_from_residual(
        graph, active_detectors, config, graph_cache, level);

    GrowingOnlyClustererResult best_result;
    ProcessingClusterBounds best_bounds{};
    AdaptiveMaxActiveCandidateScore best_max_active_score;
    bool have_best = false;
    size_t trials = 0;
    for (auto b : b_values) {
        cumulative_time_int d_upper = adaptive_max_diameter_for_buffer(config, level, b);
        if (d_upper < d_lower || d_upper <= 0) {
            continue;
        }
        ProcessingClusterBounds upper_bounds{d_upper, b};
        adaptive_fill_certificate_metadata(config, level, upper_bounds);
        if (!adaptive_bounds_satisfy_optional_constraints(config, level, upper_bounds)) {
            continue;
        }
        auto upper_config = adaptive_trial_config_with_bounds(config, level, upper_bounds);
        auto upper_result = run_growing_only_clustering_for_adaptive_trial(graph, active_detectors, level, upper_config, graph_cache);
        trials++;
        if (upper_result.accepted_clusters.empty() &&
            !(config.adaptive_defer_heavy_clusters &&
              upper_result.stats.adaptive_deferred_clusters > 0 &&
              lookahead_depth_override > 0)) {
            continue;
        }

        cumulative_time_int required_d = adaptive_required_diameter_for_result(upper_result, d_lower);
        if (required_d > d_upper) {
            continue;
        }
        ProcessingClusterBounds exact_bounds{required_d, b};
        adaptive_fill_certificate_metadata(config, level, exact_bounds);
        if (!adaptive_bounds_satisfy_optional_constraints(config, level, exact_bounds)) {
            continue;
        }
        auto exact_config = adaptive_trial_config_with_bounds(config, level, exact_bounds);
        auto exact_result = run_growing_only_clustering_for_adaptive_trial(graph, active_detectors, level, exact_config, graph_cache);
        bool take_candidate = false;
        if (config.adaptive_objective_minimize_max_active || config.adaptive_objective_minimize_max_events || config.adaptive_objective_minimize_ideal_events) {
            ProcessingClusterConfig score_config = config;
            ensure_adaptive_explicit_bounds_vector(score_config);
            ProcessingClusterBounds score_bounds = exact_bounds;
            adaptive_fill_certificate_metadata(score_config, level, score_bounds);
            score_config.explicit_bounds_by_level[level - 1] = score_bounds;
            auto score = adaptive_max_active_score_for_candidate(
                graph,
                exact_result,
                score_bounds,
                level,
                score_config,
                graph_cache,
                lookahead_depth_override);
            take_candidate = adaptive_max_active_score_is_better(score, best_max_active_score, have_best, config.adaptive_objective_minimize_max_events, config.adaptive_objective_minimize_ideal_events);
            if (take_candidate) {
                best_max_active_score = score;
            }
        } else {
            take_candidate = adaptive_result_score_is_better(
                exact_result, exact_bounds, best_result, best_bounds, have_best);
        }
        if (take_candidate) {
            best_result = std::move(exact_result);
            best_bounds = exact_bounds;
            have_best = true;
        }
    }

    if (!have_best) {
        best_bounds = config.explicit_bounds_by_level[level - 1];
        adaptive_raise_buffer_to_optional_lower_bounds(config, level, best_bounds);
        adaptive_fill_certificate_metadata(config, level, best_bounds);
        auto fallback_config = adaptive_trial_config_with_bounds(config, level, best_bounds);
        best_result = run_growing_only_clustering_for_adaptive_trial(graph, active_detectors, level, fallback_config, graph_cache);
        if (config.adaptive_objective_minimize_max_active || config.adaptive_objective_minimize_max_events || config.adaptive_objective_minimize_ideal_events) {
            ProcessingClusterConfig score_config = config;
            ensure_adaptive_explicit_bounds_vector(score_config);
            ProcessingClusterBounds score_bounds = best_bounds;
            adaptive_fill_certificate_metadata(score_config, level, score_bounds);
            score_config.explicit_bounds_by_level[level - 1] = score_bounds;
            best_max_active_score = adaptive_max_active_score_for_candidate(
                graph,
                best_result,
                score_bounds,
                level,
                score_config,
                graph_cache,
                lookahead_depth_override);
        }
    }

    config.explicit_bounds_by_level[level - 1] = best_bounds;
    best_result.bounds = best_bounds;
    if (config.adaptive_trace_bounds) {
        std::cerr << "ADAPTIVE_CLUSTER_BOUNDS_SELECTED"
                  << " strategy=breakpoint"
                  << " level=" << level
                  << " input_active=" << active_detectors.size()
                  << " d=" << best_bounds.diameter_bound
                  << " b=" << best_bounds.buffer_bound
                  << " accepted_clusters=" << best_result.accepted_clusters.size()
                  << " accepted_active=" << accepted_active_count_for_adaptive_result(best_result)
                  << " residual_active=" << best_result.residual_active_detectors.size()
                  << " trials=" << trials
                  << " breakpoints=" << b_values.size()
                  << " enforce_monotone=" << (config.adaptive_enforce_monotone_bounds ? 1 : 0)
                  << " enforce_legacy_gap=" << (config.adaptive_enforce_legacy_gap_bound ? 1 : 0)
                  << " enforce_phi=" << (config.adaptive_enforce_phi_buffer_certificate ? 1 : 0)
                  << " min_buffer_ratio=" << config.adaptive_min_buffer_to_diameter_ratio
                  << " objective=" << adaptive_objective_label(config)
                  << " lookahead=" << lookahead_depth_override
                  << " max_active=" << max_accepted_cluster_active_for_adaptive_result(best_result)
                  << " max_events=" << max_accepted_cluster_events_for_adaptive_result(best_result, config)
                  << " deferred_clusters=" << best_result.stats.adaptive_deferred_clusters
                  << " deferred_active=" << best_result.stats.adaptive_deferred_active_detectors
                  << " defer_max_active=" << config.adaptive_defer_max_active
                  << " defer_max_events=" << config.adaptive_defer_max_events
                  << " score_worst_active=" << ((config.adaptive_objective_minimize_max_active || config.adaptive_objective_minimize_max_events || config.adaptive_objective_minimize_ideal_events) ? best_max_active_score.worst_active : 0)
                  << " score_worst_events=" << ((config.adaptive_objective_minimize_max_events || config.adaptive_objective_minimize_ideal_events) ? best_max_active_score.worst_events : 0)
                  << " score_ideal_events=" << (config.adaptive_objective_minimize_ideal_events ? best_max_active_score.ideal_events : 0)
                  << " score_weighted_ideal_events=" << (config.adaptive_objective_minimize_ideal_events ? best_max_active_score.weighted_ideal_events : 0.0)
                  << " score_nonempty_event_levels=" << (config.adaptive_objective_minimize_ideal_events ? best_max_active_score.nonempty_event_levels : 0)
                  << " score_residual_after_lookahead=" << ((config.adaptive_objective_minimize_max_active || config.adaptive_objective_minimize_max_events || config.adaptive_objective_minimize_ideal_events) ? best_max_active_score.residual_after_lookahead : 0)
                  << "\n";
    }
    return best_result;
}


GrowingOnlyClustererResult choose_and_run_adaptive_growing_only_clustering(
    const MatchingGraph& graph,
    std::vector<uint64_t> active_detectors,
    size_t level,
    ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache) {
    if (config.adaptive_breakpoint_strategy) {
        return choose_and_run_breakpoint_adaptive_growing_only_clustering(
            graph,
            std::move(active_detectors),
            level,
            config,
            graph_cache,
            config.adaptive_lookahead_levels);
    }

    ensure_adaptive_explicit_bounds_vector(config);
    std::sort(active_detectors.begin(), active_detectors.end());
    active_detectors.erase(std::unique(active_detectors.begin(), active_detectors.end()), active_detectors.end());

    const auto d_values = adaptive_candidate_values_from_residual(
        graph, active_detectors, config, graph_cache, level, false);
    const auto b_values = adaptive_candidate_values_from_residual(
        graph, active_detectors, config, graph_cache, level, true);

    GrowingOnlyClustererResult best_result;
    ProcessingClusterBounds best_bounds{};
    bool have_best = false;
    size_t trials = 0;
    const size_t max_trials = std::max<size_t>(1, config.adaptive_candidate_limit);
    for (auto d : d_values) {
        for (auto b : b_values) {
            if (trials >= max_trials) {
                break;
            }
            ProcessingClusterBounds candidate_bounds{d, b};
            if (!adaptive_raise_buffer_to_optional_lower_bounds(config, level, candidate_bounds) ||
                !adaptive_bounds_satisfy_optional_constraints(config, level, candidate_bounds)) {
                continue;
            }
            trials++;
            auto trial_config = adaptive_trial_config_with_bounds(config, level, candidate_bounds);
            auto trial_result = run_growing_only_clustering_for_adaptive_trial(
                graph, active_detectors, level, trial_config, graph_cache);
            if (adaptive_result_score_is_better(
                    trial_result, candidate_bounds, best_result, best_bounds, have_best)) {
                best_result = std::move(trial_result);
                best_bounds = candidate_bounds;
                have_best = true;
            }
        }
        if (trials >= max_trials) {
            break;
        }
    }

    if (!have_best) {
        // Fall back to the current scheduled level bound if every adaptive
        // candidate was removed by optional d/b constraints.
        best_bounds = config.explicit_bounds_by_level[level - 1];
        auto fallback_config = adaptive_trial_config_with_bounds(config, level, best_bounds);
        best_result = run_growing_only_clustering_for_adaptive_trial(graph, active_detectors, level, fallback_config, graph_cache);
    }

    config.explicit_bounds_by_level[level - 1] = best_bounds;
    best_result.bounds = best_bounds;
    if (config.adaptive_trace_bounds) {
        std::cerr << "ADAPTIVE_CLUSTER_BOUNDS_SELECTED"
                  << " level=" << level
                  << " input_active=" << active_detectors.size()
                  << " d=" << best_bounds.diameter_bound
                  << " b=" << best_bounds.buffer_bound
                  << " accepted_clusters=" << best_result.accepted_clusters.size()
                  << " accepted_active=" << accepted_active_count_for_adaptive_result(best_result)
                  << " residual_active=" << best_result.residual_active_detectors.size()
                  << " trials=" << trials
                  << " enforce_monotone=" << (config.adaptive_enforce_monotone_bounds ? 1 : 0)
                  << " enforce_legacy_gap=" << (config.adaptive_enforce_legacy_gap_bound ? 1 : 0)
                  << " enforce_phi=" << (config.adaptive_enforce_phi_buffer_certificate ? 1 : 0)
                  << " min_buffer_ratio=" << config.adaptive_min_buffer_to_diameter_ratio
                  << "\n";
    }
    return best_result;
}

GrowingOnlyClustererResult run_growing_only_clustering(
    const MatchingGraph& graph,
    std::vector<uint64_t> active_detectors,
    size_t level,
    const ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache) {
    GrowingOnlyClusterer clusterer(graph, std::move(active_detectors), level, config, graph_cache);
    return clusterer.run();
}


GrowingOnlyClustererResult run_growing_only_clustering_with_adaptive_bounds(
    const MatchingGraph& graph,
    std::vector<uint64_t> active_detectors,
    size_t level,
    ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache) {
    if (!config.adaptive_cluster_bounds) {
        return run_growing_only_clustering(graph, std::move(active_detectors), level, config, graph_cache);
    }
    return choose_and_run_adaptive_growing_only_clustering(
        graph, std::move(active_detectors), level, config, graph_cache);
}

}  // namespace pm
