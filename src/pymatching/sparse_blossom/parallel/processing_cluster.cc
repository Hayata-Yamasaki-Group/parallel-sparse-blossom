// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "pymatching/sparse_blossom/parallel/processing_cluster.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <future>
#include <thread>
#include <chrono>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
#include <stdexcept>
#include <sstream>
#include <string>
#include <unordered_set>

namespace pm {
namespace {

using steady_clock = std::chrono::steady_clock;

uint64_t elapsed_ns(steady_clock::time_point start) {
    const char* disable = std::getenv("PYMATCHING_DISABLE_WALL_TIMING");
    const char* event_only_disable = std::getenv("PYMATCHING_EVENTCOUNT_ONLY_NO_WALL_TIMING");
    if ((disable != nullptr && std::string(disable) == "1") ||
        (event_only_disable != nullptr && std::string(event_only_disable) == "1")) {
        return 0;
    }
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(steady_clock::now() - start).count());
}

constexpr cumulative_time_int INF_DISTANCE = std::numeric_limits<cumulative_time_int>::max() / 4;
constexpr uint64_t DEFAULT_PROCESSING_CLUSTER_CACHE_MAX_BYTES = 2048ull * 1024ull * 1024ull;
constexpr size_t MAX_PROCESSING_CLUSTER_ALL_PAIRS_NODES = 7000;
constexpr uint32_t PACKED_INF_DISTANCE = std::numeric_limits<uint32_t>::max();

size_t processing_cluster_all_pairs_max_nodes() {
    if (const char* e = std::getenv("PYMATCHING_PROCESSING_CLUSTER_ALL_PAIRS_MAX_NODES")) {
        return static_cast<size_t>(std::stoull(e));
    }
    return MAX_PROCESSING_CLUSTER_ALL_PAIRS_NODES;
}

uint64_t processing_cluster_cache_max_bytes() {
    if (const char* e = std::getenv("PYMATCHING_PROCESSING_CLUSTER_CACHE_MAX_BYTES")) {
        return static_cast<uint64_t>(std::stoull(e));
    }
    if (const char* e = std::getenv("PYMATCHING_PROCESSING_CLUSTER_CACHE_MAX_GIB")) {
        long double gib = std::stold(e);
        if (!(gib > 0)) {
            return 0;
        }
        long double bytes = gib * 1024.0L * 1024.0L * 1024.0L;
        long double max_u64 = static_cast<long double>(std::numeric_limits<uint64_t>::max());
        if (bytes >= max_u64) {
            return std::numeric_limits<uint64_t>::max();
        }
        return static_cast<uint64_t>(bytes);
    }
    return DEFAULT_PROCESSING_CLUSTER_CACHE_MAX_BYTES;
}

cumulative_time_int clamp_to_distance(double value, bool round_up);
double level_exponent(size_t x);

bool bounds_equal(const ProcessingClusterBounds& a, const ProcessingClusterBounds& b) {
    return a.diameter_bound == b.diameter_bound && a.buffer_bound == b.buffer_bound;
}

cumulative_time_int ceil_long_double_to_distance(long double value) {
    if (!(value > 0)) {
        return 0;
    }
    auto max_distance = std::numeric_limits<cumulative_time_int>::max() / 8;
    if (value >= static_cast<long double>(max_distance)) {
        return max_distance;
    }
    long double rounded = std::round(value);
    long double tolerance = std::max(1e-9L, std::fabs(value) * 1e-15L);
    if (std::fabs(value - rounded) <= tolerance) {
        value = rounded;
    }
    return static_cast<cumulative_time_int>(std::ceil(value));
}

cumulative_time_int strict_threshold_to_distance_requirement(
    long double threshold,
    cumulative_time_int min_buffer_slack) {
    if (min_buffer_slack < 1) {
        min_buffer_slack = 1;
    }
    auto max_distance = std::numeric_limits<cumulative_time_int>::max() / 8;
    if (threshold >= static_cast<long double>(max_distance)) {
        return max_distance;
    }
    return static_cast<cumulative_time_int>(std::floor(threshold)) + min_buffer_slack;
}

ProcessingClusterBounds make_phi_certified_processing_cluster_bounds(
    cumulative_time_int diameter_bound,
    cumulative_time_int buffer_bound,
    long double phi,
    cumulative_time_int phi_required_buffer_bound,
    cumulative_time_int ratio_required_buffer_bound,
    cumulative_time_int min_required_buffer_bound,
    cumulative_time_int phi_budget_required_buffer_bound = 0,
    long double schedule_diameter_bound = -1,
    long double schedule_buffer_bound = -1) {
    ProcessingClusterBounds bounds;
    bounds.diameter_bound = diameter_bound;
    bounds.buffer_bound = buffer_bound;
    if (schedule_diameter_bound > 0 && schedule_buffer_bound > 0 &&
        std::isfinite(schedule_diameter_bound) && std::isfinite(schedule_buffer_bound)) {
        bounds.has_long_double_schedule_bounds = true;
        bounds.schedule_diameter_bound = schedule_diameter_bound;
        bounds.schedule_buffer_bound = schedule_buffer_bound;
        bounds.runtime_bounds_clamped =
            static_cast<long double>(diameter_bound) + 0.5L < schedule_diameter_bound ||
            static_cast<long double>(buffer_bound) + 0.5L < schedule_buffer_bound;
    }
    bounds.has_phi_buffer_certificate = std::isfinite(phi) && phi > 0;
    bounds.phi = bounds.has_phi_buffer_certificate ? static_cast<double>(phi) : 0;
    bounds.phi_required_buffer_bound = phi_required_buffer_bound;
    bounds.phi_budget_required_buffer_bound = phi_budget_required_buffer_bound;
    bounds.ratio_required_buffer_bound = ratio_required_buffer_bound;
    bounds.min_required_buffer_bound = min_required_buffer_bound;
    return bounds;
}

long double certificate_diameter_bound(const ProcessingClusterBounds& bounds) {
    return bounds.has_long_double_schedule_bounds
        ? bounds.schedule_diameter_bound
        : static_cast<long double>(bounds.diameter_bound);
}

long double certificate_buffer_bound(const ProcessingClusterBounds& bounds) {
    return bounds.has_long_double_schedule_bounds
        ? bounds.schedule_buffer_bound
        : static_cast<long double>(bounds.buffer_bound);
}

size_t effective_max_level_unvalidated(const ProcessingClusterConfig& config) {
    if (config.max_level == 0) {
        return 0;
    }
    if (!config.explicit_bounds_by_level.empty()) {
        return std::min(config.max_level, config.explicit_bounds_by_level.size());
    }
    return config.max_level;
}


size_t max_active_detectors_for_level_unvalidated(size_t level, const ProcessingClusterConfig& config) {
    if (level == 0) {
        throw std::invalid_argument("Processing cluster levels are one-based; level 0 is invalid.");
    }
    if (level > config.max_active_detectors_by_level.size()) {
        return 0;
    }
    return config.max_active_detectors_by_level[level - 1];
}

ProcessingClusterBounds processing_cluster_bounds_for_level_unvalidated(
    size_t level, const ProcessingClusterConfig& config) {
    if (level == 0) {
        throw std::invalid_argument("Processing cluster levels are one-based; level 0 is invalid.");
    }
    if (!config.explicit_bounds_by_level.empty()) {
        if (level > config.explicit_bounds_by_level.size()) {
            throw std::invalid_argument(
                "ProcessingClusterConfig.explicit_bounds_by_level does not contain a bound for level " +
                std::to_string(level) + ".");
        }
        return config.explicit_bounds_by_level[level - 1];
    }

    auto buffer_exponent = level_exponent(level + 1);
    auto diameter_exponent = level_exponent(level);
    auto buffer = config.beta * std::pow(config.lambda, buffer_exponent) + 1;
    auto diameter = config.gamma * std::pow(config.lambda, diameter_exponent) - 1;
    return ProcessingClusterBounds{
        clamp_to_distance(diameter, false),
        clamp_to_distance(buffer, true),
    };
}

std::vector<ProcessingClusterBounds> processing_cluster_bounds_vector_unvalidated(
    const ProcessingClusterConfig& config, size_t max_level) {
    std::vector<ProcessingClusterBounds> bounds(max_level + 1);
    for (size_t level = 1; level <= max_level; level++) {
        bounds[level] = processing_cluster_bounds_for_level_unvalidated(level, config);
    }
    return bounds;
}

bool cached_bounds_match(
    const std::vector<ProcessingClusterBounds>& cached_bounds_by_level,
    const ProcessingClusterConfig& config,
    size_t max_level) {
    if (cached_bounds_by_level.size() <= max_level) {
        return false;
    }
    for (size_t level = 1; level <= max_level; level++) {
        if (!bounds_equal(cached_bounds_by_level[level],
                          processing_cluster_bounds_for_level_unvalidated(level, config))) {
            return false;
        }
    }
    return true;
}

void validate_config(const ProcessingClusterConfig& config) {
    if (config.max_level == 0) {
        throw std::invalid_argument("ProcessingClusterConfig.max_level must be at least 1.");
    }
    if (config.explicit_bounds_by_level.empty()) {
        if (!(config.beta >= 0) || !std::isfinite(config.beta)) {
            throw std::invalid_argument("ProcessingClusterConfig.beta must be a finite non-negative value.");
        }
        if (!(config.gamma >= 0) || !std::isfinite(config.gamma)) {
            throw std::invalid_argument("ProcessingClusterConfig.gamma must be a finite non-negative value.");
        }
        if (!(config.lambda > 1) || !std::isfinite(config.lambda)) {
            throw std::invalid_argument("ProcessingClusterConfig.lambda must be a finite value greater than 1.");
        }
    } else {
        for (size_t k = 0; k < config.explicit_bounds_by_level.size(); k++) {
            const auto& bounds = config.explicit_bounds_by_level[k];
            if (bounds.diameter_bound < 0 || bounds.buffer_bound < 0) {
                throw std::invalid_argument(
                    "ProcessingClusterConfig.explicit_bounds_by_level contains a negative bound at level " +
                    std::to_string(k + 1) + ".");
            }
        }
    }
    for (size_t level = 1; level <= config.max_active_detectors_by_level.size(); level++) {
        size_t cap = config.max_active_detectors_by_level[level - 1];
        if (cap == 1) {
            throw std::invalid_argument(
                "ProcessingClusterConfig.max_active_detectors_by_level level " +
                std::to_string(level) +
                " is 1. Use 0 for no cap, or at least 2 for a real cluster-size cap.");
        }
    }
    if (effective_max_level_unvalidated(config) == 0) {
        throw std::invalid_argument(
            "ProcessingClusterConfig has no usable clustering levels.  Increase max_level or provide explicit bounds.");
    }
    if (config.enforce_stopping_lemma_bounds) {
        validate_processing_cluster_stopping_lemma_bounds(config);
    }
}

cumulative_time_int clamp_to_distance(double value, bool round_up) {
    if (!(value > 0)) {
        return 0;
    }
    auto limit = static_cast<double>(std::numeric_limits<cumulative_time_int>::max() / 8);
    if (!std::isfinite(value) || value >= limit) {
        return std::numeric_limits<cumulative_time_int>::max() / 8;
    }
    return static_cast<cumulative_time_int>(round_up ? std::ceil(value) : std::floor(value));
}

double level_exponent(size_t x) {
    if (x == 0) {
        return 0;
    }
    return static_cast<double>(x) * std::log(static_cast<double>(x));
}

std::vector<uint64_t> normalise_active_detectors(
    const MatchingGraph& graph, const std::vector<uint64_t>& active_detectors) {
    bool already_sorted_unique = true;
    bool has_user_boundary_hits = false;
    uint64_t previous = 0;
    for (size_t k = 0; k < active_detectors.size(); k++) {
        uint64_t detector = active_detectors[k];
        if (detector >= graph.nodes.size()) {
            throw std::invalid_argument(
                "The detection event with index " + std::to_string(detector) +
                " does not correspond to a node in the graph, which only has " +
                std::to_string(graph.nodes.size()) + " nodes.");
        }
        if (k != 0 && detector <= previous) {
            already_sorted_unique = false;
        }
        previous = detector;
        if (detector < graph.is_user_graph_boundary_node.size() && graph.is_user_graph_boundary_node[detector]) {
            has_user_boundary_hits = true;
        }
    }

    if (already_sorted_unique && !has_user_boundary_hits) {
        return active_detectors;
    }

    std::vector<uint64_t> result = active_detectors;
    if (!already_sorted_unique) {
        std::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
    }
    if (!has_user_boundary_hits) {
        return result;
    }

    std::vector<uint64_t> filtered;
    filtered.reserve(result.size());
    for (auto detector : result) {
        if (detector < graph.is_user_graph_boundary_node.size() && graph.is_user_graph_boundary_node[detector]) {
            continue;
        }
        filtered.push_back(detector);
    }
    return filtered;
}

std::vector<cumulative_time_int> dijkstra_from_sources(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& sources,
    cumulative_time_int max_distance = INF_DISTANCE) {
    std::vector<cumulative_time_int> distances(graph.nodes.size(), INF_DISTANCE);
    using QueueEntry = std::pair<cumulative_time_int, size_t>;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;

    for (auto source : sources) {
        distances[source] = 0;
        queue.push({0, source});
    }

    while (!queue.empty()) {
        auto [distance, node] = queue.top();
        queue.pop();
        if (distance != distances[node]) {
            continue;
        }
        if (distance > max_distance) {
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
            if (distance > INF_DISTANCE - edge_weight) {
                continue;
            }
            auto candidate = distance + edge_weight;
            if (candidate > max_distance) {
                continue;
            }
            if (candidate < distances[neighbor_index]) {
                distances[neighbor_index] = candidate;
                queue.push({candidate, neighbor_index});
            }
        }
    }

    return distances;
}

std::vector<cumulative_time_int> dijkstra_from_boundary_edges(const MatchingGraph& graph) {
    std::vector<cumulative_time_int> distances(graph.nodes.size(), INF_DISTANCE);
    using QueueEntry = std::pair<cumulative_time_int, size_t>;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;

    for (size_t node = 0; node < graph.nodes.size(); node++) {
        cumulative_time_int best_boundary_weight = INF_DISTANCE;
        for (size_t k = 0; k < graph.nodes[node].neighbors.size(); k++) {
            if (graph.nodes[node].neighbors[k] != nullptr) {
                continue;
            }
            best_boundary_weight = std::min(
                best_boundary_weight,
                static_cast<cumulative_time_int>(graph.nodes[node].neighbor_weights[k]));
        }
        if (best_boundary_weight < INF_DISTANCE) {
            distances[node] = best_boundary_weight;
            queue.push({best_boundary_weight, node});
        }
    }

    while (!queue.empty()) {
        auto [distance, node] = queue.top();
        queue.pop();
        if (distance != distances[node]) {
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
            if (distance > INF_DISTANCE - edge_weight) {
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

struct FloodedThresholdComponents {
    std::vector<std::vector<size_t>> components;
};

inline cumulative_time_int cached_distance(
    const ProcessingClusterGraphCache& graph_cache, size_t source, size_t destination);

struct DisjointSet {
    std::vector<size_t> parent;
    std::vector<uint8_t> rank;

    explicit DisjointSet(size_t n) : parent(n), rank(n, 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    size_t find(size_t k) {
        while (parent[k] != k) {
            parent[k] = parent[parent[k]];
            k = parent[k];
        }
        return k;
    }

    bool unite(size_t a, size_t b) {
        a = find(a);
        b = find(b);
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
    }
};

FloodedThresholdComponents collect_components_from_disjoint_set(DisjointSet& dsu, size_t num_items) {
    FloodedThresholdComponents result;
    if (num_items == 0) {
        return result;
    }
    std::vector<size_t> root_to_component(num_items, SIZE_MAX);
    result.components.reserve(num_items);
    for (size_t item = 0; item < num_items; item++) {
        size_t root = dsu.find(item);
        size_t component_index = root_to_component[root];
        if (component_index == SIZE_MAX) {
            component_index = result.components.size();
            root_to_component[root] = component_index;
            result.components.emplace_back();
        }
        result.components[component_index].push_back(item);
    }
    return result;
}

FloodedThresholdComponents build_threshold_components_with_lightweight_cluster_flooder(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& active_detectors,
    cumulative_time_int buffer_bound) {
    FloodedThresholdComponents result;
    if (active_detectors.empty()) {
        return result;
    }

    std::vector<cumulative_time_int> distances(graph.nodes.size(), INF_DISTANCE);
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
        if (distance > buffer_bound) {
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
            if (distance > INF_DISTANCE - edge_weight) {
                continue;
            }
            auto candidate = distance + edge_weight;
            if (candidate > buffer_bound) {
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

    DisjointSet dsu(active_detectors.size());
    for (size_t node = 0; node < graph.nodes.size(); node++) {
        auto owning_source = owner[node];
        if (owning_source == SIZE_MAX || distances[node] > buffer_bound) {
            continue;
        }
        auto self_source = source_detector_index[node];
        if (self_source >= 0 && static_cast<size_t>(self_source) != owning_source) {
            dsu.unite(owning_source, static_cast<size_t>(self_source));
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
            if (neighbor_owner == SIZE_MAX || neighbor_owner == owning_source || distances[neighbor_index] > buffer_bound) {
                continue;
            }
            auto edge_weight = static_cast<cumulative_time_int>(detector_node.neighbor_weights[k]);
            if (distances[node] > INF_DISTANCE - edge_weight ||
                distances[node] + edge_weight > INF_DISTANCE - distances[neighbor_index]) {
                continue;
            }
            auto witness_distance = distances[node] + edge_weight + distances[neighbor_index];
            if (witness_distance <= buffer_bound) {
                dsu.unite(owning_source, neighbor_owner);
            }
        }
    }

    std::vector<std::vector<size_t>> components_by_root(active_detectors.size());
    for (size_t detector_index = 0; detector_index < active_detectors.size(); detector_index++) {
        components_by_root[dsu.find(detector_index)].push_back(detector_index);
    }
    for (auto& component : components_by_root) {
        if (component.empty()) {
            continue;
        }
        std::sort(component.begin(), component.end());
        result.components.push_back(std::move(component));
    }
    return result;
}

inline cumulative_time_int cached_distance(
    const ProcessingClusterGraphCache& graph_cache, size_t source, size_t destination) {
    if (!graph_cache.has_all_pairs_distances) {
        return INF_DISTANCE;
    }
    if (source > destination) {
        std::swap(source, destination);
    }
    uint64_t index = graph_cache.packed_all_pairs_row_offsets[source] + (destination - source);
    auto packed = graph_cache.packed_all_pairs_interior_distances[index];
    return packed == PACKED_INF_DISTANCE ? INF_DISTANCE : static_cast<cumulative_time_int>(packed);
}


struct ActiveDetectorDistanceCache {
    size_t num_active = 0;
    std::vector<size_t> global_to_active_index;
    std::vector<uint32_t> distances;
    std::vector<uint64_t> row_offsets;
    bool triangular = false;

    inline cumulative_time_int distance_by_active_index(size_t a, size_t b) const {
        uint64_t index;
        if (triangular) {
            if (a > b) {
                std::swap(a, b);
            }
            index = row_offsets[a] + (b - a);
        } else {
            index = a * num_active + b;
        }
        auto packed = distances[index];
        return packed == PACKED_INF_DISTANCE ? INF_DISTANCE : static_cast<cumulative_time_int>(packed);
    }

    inline cumulative_time_int distance_by_global_detector(uint64_t detector_a, uint64_t detector_b) const {
        if (detector_a >= global_to_active_index.size() || detector_b >= global_to_active_index.size()) {
            return INF_DISTANCE;
        }
        size_t a = global_to_active_index[detector_a];
        size_t b = global_to_active_index[detector_b];
        if (a == SIZE_MAX || b == SIZE_MAX) {
            return INF_DISTANCE;
        }
        return distance_by_active_index(a, b);
    }
};

template <typename Func>
void parallel_for_rows(size_t count, Func func) {
    if (count == 0) {
        return;
    }
    size_t workers = std::min<size_t>(8, std::max<size_t>(1, std::thread::hardware_concurrency()));
    if (count < 100000 || workers <= 1) {
        for (size_t i = 0; i < count; i++) {
            func(i);
        }
        return;
    }
    workers = std::min(workers, count);
    std::vector<std::future<void>> futures;
    futures.reserve(workers);
    size_t chunk = (count + workers - 1) / workers;
    for (size_t w = 0; w < workers; w++) {
        size_t start = w * chunk;
        if (start >= count) {
            break;
        }
        size_t end = std::min(count, start + chunk);
        futures.push_back(std::async(std::launch::async, [&, start, end]() {
            for (size_t i = start; i < end; i++) {
                func(i);
            }
        }));
    }
    for (auto& f : futures) {
        f.get();
    }
}

ActiveDetectorDistanceCache build_active_detector_distance_cache(
    const ProcessingClusterGraphCache& graph_cache,
    const std::vector<uint64_t>& active_detectors) {
    ActiveDetectorDistanceCache cache;
    cache.num_active = active_detectors.size();
    cache.global_to_active_index.assign(graph_cache.num_nodes, SIZE_MAX);
    for (size_t k = 0; k < active_detectors.size(); k++) {
        if (active_detectors[k] < graph_cache.num_nodes) {
            cache.global_to_active_index[active_detectors[k]] = k;
        }
    }
    cache.triangular = cache.num_active > 512;
    if (cache.triangular) {
        uint64_t num_entries =
            (static_cast<uint64_t>(cache.num_active) * static_cast<uint64_t>(cache.num_active + 1)) / 2;
        cache.distances.assign(static_cast<size_t>(num_entries), PACKED_INF_DISTANCE);
        cache.row_offsets.resize(cache.num_active);
        uint64_t next_offset = 0;
        for (size_t i = 0; i < cache.num_active; i++) {
            cache.row_offsets[i] = next_offset;
            next_offset += static_cast<uint64_t>(cache.num_active - i);
        }
        parallel_for_rows(active_detectors.size(), [&](size_t i) {
            cache.distances[cache.row_offsets[i]] = 0;
            auto source = static_cast<size_t>(active_detectors[i]);
            for (size_t j = i + 1; j < active_detectors.size(); j++) {
                auto distance = cached_distance(graph_cache, source, static_cast<size_t>(active_detectors[j]));
                cache.distances[cache.row_offsets[i] + (j - i)] =
                    distance >= INF_DISTANCE ? PACKED_INF_DISTANCE : static_cast<uint32_t>(distance);
            }
        });
    } else {
        cache.distances.assign(cache.num_active * cache.num_active, PACKED_INF_DISTANCE);
        parallel_for_rows(active_detectors.size(), [&](size_t i) {
            cache.distances[i * cache.num_active + i] = 0;
            auto source = static_cast<size_t>(active_detectors[i]);
            for (size_t j = i + 1; j < active_detectors.size(); j++) {
                auto destination = static_cast<size_t>(active_detectors[j]);
                auto distance = cached_distance(graph_cache, source, destination);
                auto packed = distance >= INF_DISTANCE ? PACKED_INF_DISTANCE : static_cast<uint32_t>(distance);
                cache.distances[i * cache.num_active + j] = packed;
                cache.distances[j * cache.num_active + i] = packed;
            }
        });
    }
    return cache;
}

FloodedThresholdComponents build_threshold_components_from_active_distance_cache(
    const ActiveDetectorDistanceCache& active_cache,
    const std::vector<uint64_t>& residual_detectors,
    cumulative_time_int buffer_bound) {
    FloodedThresholdComponents result;
    if (residual_detectors.empty()) {
        return result;
    }
    DisjointSet dsu(residual_detectors.size());
    std::vector<size_t> residual_active_indices(residual_detectors.size(), SIZE_MAX);
    for (size_t k = 0; k < residual_detectors.size(); k++) {
        auto detector = residual_detectors[k];
        if (detector < active_cache.global_to_active_index.size()) {
            residual_active_indices[k] = active_cache.global_to_active_index[detector];
        }
    }

    const uint32_t packed_buffer = buffer_bound >= static_cast<cumulative_time_int>(PACKED_INF_DISTANCE)
                                       ? PACKED_INF_DISTANCE - 1
                                       : static_cast<uint32_t>(buffer_bound);
    size_t components_remaining = residual_detectors.size();
    if (active_cache.triangular) {
        // residual_detectors is sorted and active_cache was built from the sorted
        // initial residual set, so residual_active_indices is monotone. Avoid the
        // per-pair branch/swap in distance_by_active_index on the hot O(m^2) path.
        for (size_t i = 0; i < residual_active_indices.size(); i++) {
            auto ai = residual_active_indices[i];
            if (ai == SIZE_MAX) {
                continue;
            }
            const uint64_t row_offset = active_cache.row_offsets[ai];
            for (size_t j = i + 1; j < residual_active_indices.size(); j++) {
                auto aj = residual_active_indices[j];
                if (aj == SIZE_MAX || aj < ai) {
                    continue;
                }
                auto packed = active_cache.distances[row_offset + (aj - ai)];
                if (packed <= packed_buffer && dsu.unite(i, j)) {
                    if (--components_remaining == 1) {
                        break;
                    }
                }
            }
            if (components_remaining == 1) {
                break;
            }
        }
    } else {
        for (size_t i = 0; i < residual_active_indices.size(); i++) {
            auto ai = residual_active_indices[i];
            if (ai == SIZE_MAX) {
                continue;
            }
            const uint64_t row_offset = static_cast<uint64_t>(ai) * active_cache.num_active;
            for (size_t j = i + 1; j < residual_active_indices.size(); j++) {
                auto aj = residual_active_indices[j];
                if (aj == SIZE_MAX) {
                    continue;
                }
                auto packed = active_cache.distances[row_offset + aj];
                if (packed <= packed_buffer && dsu.unite(i, j)) {
                    if (--components_remaining == 1) {
                        break;
                    }
                }
            }
            if (components_remaining == 1) {
                break;
            }
        }
    }

    std::vector<std::vector<size_t>> components_by_root(residual_detectors.size());
    for (size_t detector_index = 0; detector_index < residual_detectors.size(); detector_index++) {
        components_by_root[dsu.find(detector_index)].push_back(detector_index);
    }
    for (auto& component : components_by_root) {
        if (component.empty()) {
            continue;
        }
        result.components.push_back(std::move(component));
    }
    return result;
}

cumulative_time_int component_diameter_from_active_distance_cache_with_cutoff(
    const ActiveDetectorDistanceCache& active_cache,
    const std::vector<uint64_t>& component_detectors,
    std::optional<cumulative_time_int> cutoff) {
    cumulative_time_int diameter = 0;
    uint32_t packed_cutoff = PACKED_INF_DISTANCE - 1;
    if (cutoff.has_value() && *cutoff < static_cast<cumulative_time_int>(PACKED_INF_DISTANCE)) {
        packed_cutoff = static_cast<uint32_t>(*cutoff);
    }

    std::vector<size_t> active_indices;
    active_indices.reserve(component_detectors.size());
    for (auto detector : component_detectors) {
        if (detector >= active_cache.global_to_active_index.size()) {
            return INF_DISTANCE;
        }
        auto active_index = active_cache.global_to_active_index[detector];
        if (active_index == SIZE_MAX) {
            return INF_DISTANCE;
        }
        active_indices.push_back(active_index);
    }

    if (active_cache.triangular) {
        for (size_t i = 0; i < active_indices.size(); i++) {
            auto ai = active_indices[i];
            for (size_t j = i + 1; j < active_indices.size(); j++) {
                auto aj = active_indices[j];
                auto lo = ai;
                auto hi = aj;
                if (lo > hi) {
                    std::swap(lo, hi);
                }
                auto packed = active_cache.distances[active_cache.row_offsets[lo] + (hi - lo)];
                if (packed == PACKED_INF_DISTANCE || packed > packed_cutoff) {
                    return INF_DISTANCE;
                }
                diameter = std::max(diameter, static_cast<cumulative_time_int>(packed));
            }
        }
    } else {
        for (size_t i = 0; i < active_indices.size(); i++) {
            auto ai = active_indices[i];
            const uint64_t row_offset = static_cast<uint64_t>(ai) * active_cache.num_active;
            for (size_t j = i + 1; j < active_indices.size(); j++) {
                auto packed = active_cache.distances[row_offset + active_indices[j]];
                if (packed == PACKED_INF_DISTANCE || packed > packed_cutoff) {
                    return INF_DISTANCE;
                }
                diameter = std::max(diameter, static_cast<cumulative_time_int>(packed));
            }
        }
    }
    return diameter;
}


FloodedThresholdComponents build_single_threshold_component(size_t num_active_detectors) {
    FloodedThresholdComponents result;
    if (num_active_detectors == 0) {
        return result;
    }
    result.components.emplace_back();
    auto& component = result.components.back();
    component.resize(num_active_detectors);
    std::iota(component.begin(), component.end(), 0);
    return result;
}

bool buffer_covers_all_interior_distances(
    const ProcessingClusterGraphCache& graph_cache, cumulative_time_int buffer_bound) {
    return graph_cache.has_all_pairs_distances &&
        graph_cache.all_interior_pairs_finite &&
        graph_cache.max_finite_interior_distance < INF_DISTANCE &&
        buffer_bound >= graph_cache.max_finite_interior_distance;
}

bool buffer_covers_all_interior_distances_by_upper_bound(
    const ProcessingClusterGraphCache& graph_cache, cumulative_time_int buffer_bound) {
    return graph_cache.interior_diameter_upper_bound < INF_DISTANCE &&
        buffer_bound >= graph_cache.interior_diameter_upper_bound;
}

bool diameter_bound_covers_all_interior_distances(
    const ProcessingClusterGraphCache& graph_cache, cumulative_time_int diameter_bound) {
    return graph_cache.has_all_pairs_distances &&
        graph_cache.all_interior_pairs_finite &&
        graph_cache.max_finite_interior_distance < INF_DISTANCE &&
        diameter_bound >= graph_cache.max_finite_interior_distance;
}

FloodedThresholdComponents build_threshold_components_from_cache(
    const ProcessingClusterGraphCache& graph_cache,
    const std::vector<uint64_t>& active_detectors,
    cumulative_time_int buffer_bound) {
    FloodedThresholdComponents result;
    if (active_detectors.empty()) {
        return result;
    }

    DisjointSet dsu(active_detectors.size());
    size_t components_remaining = active_detectors.size();
    for (size_t i = 0; i < active_detectors.size(); i++) {
        auto source = static_cast<size_t>(active_detectors[i]);
        for (size_t j = i + 1; j < active_detectors.size(); j++) {
            auto destination = static_cast<size_t>(active_detectors[j]);
            if (cached_distance(graph_cache, source, destination) <= buffer_bound && dsu.unite(i, j)) {
                if (--components_remaining == 1) {
                    break;
                }
            }
        }
        if (components_remaining == 1) {
            break;
        }
    }

    return collect_components_from_disjoint_set(dsu, active_detectors.size());
}

const ProcessingClusterRadiusNeighborList* radius_neighbors_for_level(
    const ProcessingClusterGraphCache& graph_cache,
    size_t level,
    cumulative_time_int buffer_bound) {
    if (!graph_cache.has_radius_neighbor_precompute || level == 0 ||
        level >= graph_cache.radius_neighbors_by_level.size()) {
        return nullptr;
    }
    const auto& neighbors = graph_cache.radius_neighbors_by_level[level];
    if (neighbors.offsets.empty() || neighbors.buffer_bound != buffer_bound) {
        return nullptr;
    }
    return &neighbors;
}

FloodedThresholdComponents build_threshold_components_from_radius_neighbors(
    const ProcessingClusterGraphCache& graph_cache,
    const ProcessingClusterRadiusNeighborList& radius_neighbors,
    const std::vector<uint64_t>& active_detectors) {
    FloodedThresholdComponents result;
    if (active_detectors.empty()) {
        return result;
    }

    static thread_local std::vector<uint32_t> detector_active_stamp;
    static thread_local std::vector<size_t> detector_active_index;
    static thread_local uint32_t detector_active_epoch = 1;
    static thread_local std::vector<uint8_t> active_visited;
    static thread_local std::vector<size_t> stack;

    if (detector_active_stamp.size() != graph_cache.num_nodes) {
        detector_active_stamp.assign(graph_cache.num_nodes, 0);
        detector_active_index.resize(graph_cache.num_nodes);
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
            detector_active_index[detector] = k;
        }
    }

    active_visited.assign(active_detectors.size(), 0);
    stack.clear();
    result.components.reserve(std::min<size_t>(active_detectors.size(), 32));

    for (size_t start_index = 0; start_index < active_detectors.size(); start_index++) {
        if (active_visited[start_index]) {
            continue;
        }
        result.components.emplace_back();
        auto& component = result.components.back();
        component.reserve(4);
        active_visited[start_index] = 1;
        stack.push_back(start_index);
        while (!stack.empty()) {
            size_t i = stack.back();
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
                size_t j = detector_active_index[neighbor];
                if (active_visited[j]) {
                    continue;
                }
                active_visited[j] = 1;
                stack.push_back(j);
            }
        }
    }
    return result;
}

std::vector<ProcessingClusterRadiusNeighborList> build_radius_neighbor_lists_for_levels_from_all_pairs(
    const ProcessingClusterGraphCache& graph_cache,
    const ProcessingClusterConfig& config,
    size_t max_level) {
    std::vector<cumulative_time_int> buffer_bounds(max_level + 1, 0);
    for (size_t level = 1; level <= max_level; level++) {
        buffer_bounds[level] = processing_cluster_bounds_for_level(level, config).buffer_bound;
    }

    std::vector<ProcessingClusterRadiusNeighborList> results(max_level + 1);
    for (size_t level = 1; level <= max_level; level++) {
        results[level].buffer_bound = buffer_bounds[level];
        results[level].offsets.assign(graph_cache.num_nodes + 1, 0);
    }

    // Count and fill all requested levels in two triangular all-pairs passes.
    // The old code scanned the full N x N distance matrix once per level and
    // created one vector per row.  This version reads each unordered pair once,
    // writes it to every threshold that contains it, and directly materializes
    // CSR storage.  For the canonical d=13 case this removes most of the
    // one-time radius-neighbor precompute overhead while preserving the exact
    // one-based level table used by the online clustering path.
    for (size_t source = 0; source < graph_cache.num_nodes; source++) {
        for (size_t destination = source + 1; destination < graph_cache.num_nodes; destination++) {
            auto distance = cached_distance(graph_cache, source, destination);
            if (distance >= INF_DISTANCE) {
                continue;
            }
            size_t first_level = 1;
            while (first_level <= max_level && distance > buffer_bounds[first_level]) {
                first_level++;
            }
            for (size_t level = first_level; level <= max_level; level++) {
                results[level].offsets[source + 1]++;
                results[level].offsets[destination + 1]++;
            }
        }
    }

    for (size_t level = 1; level <= max_level; level++) {
        auto& offsets = results[level].offsets;
        for (size_t node = 0; node < graph_cache.num_nodes; node++) {
            offsets[node + 1] += offsets[node];
        }
        results[level].neighbors.resize(static_cast<size_t>(offsets.back()));
    }

    std::vector<std::vector<uint64_t>> cursors(max_level + 1);
    for (size_t level = 1; level <= max_level; level++) {
        cursors[level] = results[level].offsets;
    }

    for (size_t source = 0; source < graph_cache.num_nodes; source++) {
        for (size_t destination = source + 1; destination < graph_cache.num_nodes; destination++) {
            auto distance = cached_distance(graph_cache, source, destination);
            if (distance >= INF_DISTANCE) {
                continue;
            }
            size_t first_level = 1;
            while (first_level <= max_level && distance > buffer_bounds[first_level]) {
                first_level++;
            }
            for (size_t level = first_level; level <= max_level; level++) {
                results[level].neighbors[static_cast<size_t>(cursors[level][source]++)] =
                    static_cast<uint32_t>(destination);
                results[level].neighbors[static_cast<size_t>(cursors[level][destination]++)] =
                    static_cast<uint32_t>(source);
            }
        }
    }

    return results;
}

bool node_has_boundary_edge_within_buffer(
    const DetectorNode& node, cumulative_time_int distance_to_node, cumulative_time_int buffer_bound) {
    if (distance_to_node >= INF_DISTANCE) {
        return false;
    }
    for (size_t k = 0; k < node.neighbors.size(); k++) {
        if (node.neighbors[k] == nullptr) {
            auto edge_weight = static_cast<cumulative_time_int>(node.neighbor_weights[k]);
            return distance_to_node <= buffer_bound && edge_weight <= buffer_bound - distance_to_node;
        }
    }
    return false;
}

void update_boundary_reachability_and_influence(
    const MatchingGraph& graph,
    const std::vector<cumulative_time_int>& distances,
    ProcessingCluster& cluster) {
    cluster.nearest_boundary_match_distance = INF_DISTANCE;
    for (size_t node = 0; node < graph.nodes.size(); node++) {
        if (distances[node] <= cluster.buffer_bound) {
            cluster.influence_vertices.push_back(node);
        }
        if (distances[node] >= INF_DISTANCE) {
            continue;
        }
        for (size_t k = 0; k < graph.nodes[node].neighbors.size(); k++) {
            if (graph.nodes[node].neighbors[k] != nullptr) {
                continue;
            }
            auto edge_weight = static_cast<cumulative_time_int>(graph.nodes[node].neighbor_weights[k]);
            if (distances[node] > cluster.buffer_bound || edge_weight > cluster.buffer_bound - distances[node]) {
                continue;
            }
            cluster.boundary_endpoint_vertices.push_back(node);
            cluster.nearest_boundary_match_distance =
                std::min(cluster.nearest_boundary_match_distance, distances[node] + edge_weight);
            break;
        }
    }
}

void fill_influence_region(const MatchingGraph& graph, ProcessingCluster& cluster) {
    auto distances = dijkstra_from_sources(graph, cluster.active_detectors, cluster.buffer_bound);
    update_boundary_reachability_and_influence(graph, distances, cluster);
}

void fill_influence_region_from_cache(const ProcessingClusterGraphCache& graph_cache, ProcessingCluster& cluster) {
    cluster.nearest_boundary_match_distance = INF_DISTANCE;
    for (size_t node = 0; node < graph_cache.num_nodes; node++) {
        cumulative_time_int distance = INF_DISTANCE;
        for (auto source : cluster.active_detectors) {
            distance = std::min(distance, cached_distance(graph_cache, source, node));
        }
        if (distance <= cluster.buffer_bound) {
            cluster.influence_vertices.push_back(node);
        }
        auto boundary_weight = graph_cache.boundary_edge_weight_by_vertex[node];
        if (distance > cluster.buffer_bound || boundary_weight >= INF_DISTANCE) {
            continue;
        }
        if (boundary_weight > cluster.buffer_bound - distance) {
            continue;
        }
        cluster.boundary_endpoint_vertices.push_back(node);
        cluster.nearest_boundary_match_distance =
            std::min(cluster.nearest_boundary_match_distance, distance + boundary_weight);
    }
}

void fill_influence_region_from_cache(
    const ProcessingClusterGraphCache& graph_cache,
    ProcessingCluster& cluster,
    std::vector<cumulative_time_int>& min_distances) {
    std::fill(min_distances.begin(), min_distances.end(), INF_DISTANCE);
    for (auto source : cluster.active_detectors) {
        for (size_t node = 0; node < graph_cache.num_nodes; node++) {
            min_distances[node] = std::min(min_distances[node], cached_distance(graph_cache, source, node));
        }
    }

    cluster.nearest_boundary_match_distance = INF_DISTANCE;
    cluster.influence_vertices.clear();
    cluster.boundary_endpoint_vertices.clear();
    for (size_t node = 0; node < graph_cache.num_nodes; node++) {
        auto distance = min_distances[node];
        if (distance <= cluster.buffer_bound) {
            cluster.influence_vertices.push_back(node);
        }
        auto boundary_weight = graph_cache.boundary_edge_weight_by_vertex[node];
        if (distance > cluster.buffer_bound || boundary_weight >= INF_DISTANCE) {
            continue;
        }
        if (boundary_weight > cluster.buffer_bound - distance) {
            continue;
        }
        cluster.boundary_endpoint_vertices.push_back(node);
        cluster.nearest_boundary_match_distance =
            std::min(cluster.nearest_boundary_match_distance, distance + boundary_weight);
    }
}

void fill_influence_region_from_radius_neighbors(
    const ProcessingClusterGraphCache& graph_cache,
    const ProcessingClusterRadiusNeighborList& radius_neighbors,
    ProcessingCluster& cluster) {
    static thread_local std::vector<uint32_t> influence_stamp;
    static thread_local uint32_t influence_epoch = 1;
    if (influence_stamp.size() != graph_cache.num_nodes) {
        influence_stamp.assign(graph_cache.num_nodes, 0);
        influence_epoch = 1;
    } else if (influence_epoch == std::numeric_limits<uint32_t>::max()) {
        std::fill(influence_stamp.begin(), influence_stamp.end(), 0);
        influence_epoch = 1;
    } else {
        influence_epoch++;
    }

    cluster.influence_vertices.clear();
    cluster.boundary_endpoint_vertices.clear();
    cluster.nearest_boundary_match_distance = INF_DISTANCE;

    auto mark_vertex = [&](size_t node) {
        if (node >= influence_stamp.size() || influence_stamp[node] == influence_epoch) {
            return;
        }
        influence_stamp[node] = influence_epoch;
        cluster.influence_vertices.push_back(node);
    };

    for (auto source64 : cluster.active_detectors) {
        size_t source = static_cast<size_t>(source64);
        if (source >= graph_cache.num_nodes || source + 1 >= radius_neighbors.offsets.size()) {
            continue;
        }
        mark_vertex(source);
        if (source < graph_cache.nearest_boundary_match_distance_by_vertex.size()) {
            cluster.nearest_boundary_match_distance = std::min(
                cluster.nearest_boundary_match_distance,
                graph_cache.nearest_boundary_match_distance_by_vertex[source]);
        }
        uint64_t begin = radius_neighbors.offsets[source];
        uint64_t end = radius_neighbors.offsets[source + 1];
        for (uint64_t p = begin; p < end; p++) {
            mark_vertex(radius_neighbors.neighbors[p]);
        }
    }

    // Subgraph construction does not require sorted influence vertices, and
    // nearest_boundary_match_distance has already been computed from the graph
    // cache for the active sources. Avoid sorting and duplicate boundary endpoint
    // reconstruction in this hot path.
}



void fill_full_interior_influence_region_from_cache(
    const ProcessingClusterGraphCache& graph_cache, ProcessingCluster& cluster) {
    cluster.influence_vertices.clear();
    cluster.influence_is_full_graph = true;
    // Full-graph identity subgraphs do not need an explicit influence vector or
    // boundary endpoint list online; the global detector graph is already shared.
    cluster.boundary_endpoint_vertices.clear();
    cluster.nearest_boundary_match_distance = INF_DISTANCE;
    if (cluster.active_detectors.size() % 2 == 1) {
        for (auto source : cluster.active_detectors) {
            if (source < graph_cache.nearest_boundary_match_distance_by_vertex.size()) {
                cluster.nearest_boundary_match_distance = std::min(
                    cluster.nearest_boundary_match_distance,
                    graph_cache.nearest_boundary_match_distance_by_vertex[source]);
            }
        }
    }
}

void fill_full_graph_influence_region_from_cache(
    const ProcessingClusterGraphCache& graph_cache, ProcessingCluster& cluster) {
    cluster.influence_vertices.clear();
    cluster.influence_is_full_graph = true;
    cluster.boundary_endpoint_vertices.clear();
    cluster.nearest_boundary_match_distance = INF_DISTANCE;
    if (cluster.active_detectors.size() % 2 == 1) {
        for (auto source : cluster.active_detectors) {
            if (source >= graph_cache.nearest_boundary_match_distance_by_vertex.size()) {
                continue;
            }
            cluster.nearest_boundary_match_distance = std::min(
                cluster.nearest_boundary_match_distance,
                graph_cache.nearest_boundary_match_distance_by_vertex[source]);
        }
    }
}

bool can_use_full_graph_influence_fast_path(
    const ProcessingClusterGraphCache& graph_cache, cumulative_time_int buffer_bound) {
    cumulative_time_int interior_bound = INF_DISTANCE;
    if (graph_cache.has_all_pairs_distances &&
        graph_cache.all_interior_pairs_finite &&
        graph_cache.max_finite_interior_distance < INF_DISTANCE) {
        interior_bound = graph_cache.max_finite_interior_distance;
    } else if (graph_cache.interior_diameter_upper_bound < INF_DISTANCE) {
        // Safe graph-only upper bound from landmark Dijkstras. This may be loose,
        // but if it fits inside the buffer then the full graph is certainly in
        // the influence region without any shot-dependent Dijkstra.
        interior_bound = graph_cache.interior_diameter_upper_bound;
    }
    if (interior_bound >= INF_DISTANCE || graph_cache.max_boundary_edge_weight >= INF_DISTANCE) {
        return false;
    }
    if (interior_bound > INF_DISTANCE - graph_cache.max_boundary_edge_weight) {
        return false;
    }
    return buffer_bound >= interior_bound + graph_cache.max_boundary_edge_weight;
}

cumulative_time_int component_diameter_from_cache_with_cutoff(
    const ProcessingClusterGraphCache& graph_cache,
    const std::vector<uint64_t>& component_detectors,
    std::optional<cumulative_time_int> cutoff) {
    cumulative_time_int diameter = 0;
    for (size_t i = 0; i < component_detectors.size(); i++) {
        for (size_t j = i + 1; j < component_detectors.size(); j++) {
            auto distance = cached_distance(graph_cache, component_detectors[i], component_detectors[j]);
            if (distance >= INF_DISTANCE) {
                return INF_DISTANCE;
            }
            if (cutoff.has_value() && distance > *cutoff) {
                return INF_DISTANCE;
            }
            diameter = std::max(diameter, distance);
        }
    }
    return diameter;
}


cumulative_time_int nearest_boundary_match_distance_for_sources_from_cache(
    const ProcessingClusterGraphCache& graph_cache,
    const std::vector<uint64_t>& component_detectors) {
    cumulative_time_int best = INF_DISTANCE;
    for (auto source : component_detectors) {
        if (source < graph_cache.nearest_boundary_match_distance_by_vertex.size()) {
            best = std::min(best, graph_cache.nearest_boundary_match_distance_by_vertex[source]);
        }
    }
    return best;
}

cumulative_time_int nearest_boundary_match_distance_for_sources(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& component_detectors,
    const ProcessingClusterGraphCache* graph_cache) {
    if (component_detectors.empty()) {
        return INF_DISTANCE;
    }
    if (graph_cache != nullptr) {
        return nearest_boundary_match_distance_for_sources_from_cache(*graph_cache, component_detectors);
    }
    auto distances = dijkstra_from_sources(graph, component_detectors);
    cumulative_time_int best = INF_DISTANCE;
    for (size_t node = 0; node < graph.nodes.size(); node++) {
        if (distances[node] >= INF_DISTANCE) {
            continue;
        }
        for (size_t k = 0; k < graph.nodes[node].neighbors.size(); k++) {
            if (graph.nodes[node].neighbors[k] != nullptr) {
                continue;
            }
            auto edge_weight = static_cast<cumulative_time_int>(graph.nodes[node].neighbor_weights[k]);
            if (distances[node] > INF_DISTANCE - edge_weight) {
                continue;
            }
            best = std::min(best, distances[node] + edge_weight);
        }
    }
    return best;
}

cumulative_time_int farthest_boundary_match_distance_for_sources_from_cache(
    const ProcessingClusterGraphCache& graph_cache, const std::vector<uint64_t>& component_detectors) {
    if (component_detectors.empty()) {
        return INF_DISTANCE;
    }
    cumulative_time_int worst = 0;
    for (auto source : component_detectors) {
        if (source >= graph_cache.nearest_boundary_match_distance_by_vertex.size()) {
            return INF_DISTANCE;
        }
        worst = std::max(worst, graph_cache.nearest_boundary_match_distance_by_vertex[source]);
    }
    return worst;
}

cumulative_time_int farthest_boundary_match_distance_for_sources(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& component_detectors,
    const ProcessingClusterGraphCache* graph_cache) {
    if (component_detectors.empty()) {
        return INF_DISTANCE;
    }
    if (graph_cache != nullptr) {
        return farthest_boundary_match_distance_for_sources_from_cache(*graph_cache, component_detectors);
    }
    auto boundary_distances = dijkstra_from_boundary_edges(graph);
    cumulative_time_int worst = 0;
    for (auto detector : component_detectors) {
        if (detector >= boundary_distances.size()) {
            return INF_DISTANCE;
        }
        worst = std::max(worst, boundary_distances[detector]);
    }
    return worst;
}

bool cluster_can_stop_locally_from_parity_and_boundary(
    size_t active_detector_count,
    cumulative_time_int boundary_inclusive_diameter,
    cumulative_time_int diameter_bound) {
    if (active_detector_count % 2 == 0) {
        return true;
    }
    return boundary_inclusive_diameter <= diameter_bound;
}

bool cluster_can_stop_locally(const ProcessingCluster& cluster) {
    if (cluster.active_detectors.size() % 2 == 0) {
        return true;
    }
    if (cluster.nearest_boundary_match_distance > cluster.diameter_bound) {
        return false;
    }
    return true;
}

cumulative_time_int component_diameter_within_bound(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& component_detectors,
    cumulative_time_int diameter_bound) {
    cumulative_time_int diameter = 0;
    for (auto detector : component_detectors) {
        auto from_source = dijkstra_from_sources(graph, {detector}, diameter_bound);
        for (auto other : component_detectors) {
            auto distance = from_source[other];
            if (distance > diameter_bound) {
                return INF_DISTANCE;
            }
            diameter = std::max(diameter, distance);
        }
    }

    return diameter;
}

cumulative_time_int component_diameter_from_cache(
    const ProcessingClusterGraphCache& graph_cache, const std::vector<uint64_t>& component_detectors) {
    cumulative_time_int diameter = 0;
    for (auto detector_a : component_detectors) {
        for (auto detector_b : component_detectors) {
            auto distance = cached_distance(graph_cache, detector_a, detector_b);
            if (distance >= INF_DISTANCE) {
                return INF_DISTANCE;
            }
            diameter = std::max(diameter, distance);
        }
    }
    return diameter;
}


void fill_interior_component_upper_bounds(const MatchingGraph& graph, ProcessingClusterGraphCache& graph_cache) {
    graph_cache.interior_component_id_by_vertex.assign(graph.nodes.size(), SIZE_MAX);
    graph_cache.interior_component_diameter_upper_bounds.clear();
    std::vector<size_t> stack;
    std::vector<size_t> component_nodes;
    for (size_t start = 0; start < graph.nodes.size(); start++) {
        if (graph_cache.interior_component_id_by_vertex[start] != SIZE_MAX) {
            continue;
        }
        size_t component_id = graph_cache.interior_component_diameter_upper_bounds.size();
        component_nodes.clear();
        stack.clear();
        stack.push_back(start);
        graph_cache.interior_component_id_by_vertex[start] = component_id;
        while (!stack.empty()) {
            size_t node = stack.back();
            stack.pop_back();
            component_nodes.push_back(node);
            const auto& detector_node = graph.nodes[node];
            for (auto neighbor : detector_node.neighbors) {
                if (neighbor == nullptr) {
                    continue;
                }
                size_t neighbor_id = static_cast<size_t>(neighbor - graph.nodes.data());
                if (graph_cache.interior_component_id_by_vertex[neighbor_id] != SIZE_MAX) {
                    continue;
                }
                graph_cache.interior_component_id_by_vertex[neighbor_id] = component_id;
                stack.push_back(neighbor_id);
            }
        }

        if (const char* e = std::getenv("PYMATCHING_COMPONENT_IDS_ONLY_GRAPH_CACHE")) {
            if (std::string(e) == "1") {
                graph_cache.interior_component_diameter_upper_bounds.push_back(INF_DISTANCE);
                continue;
            }
        }

        std::vector<size_t> landmark_candidates;
        landmark_candidates.push_back(component_nodes.front());
        landmark_candidates.push_back(component_nodes[component_nodes.size() / 2]);
        landmark_candidates.push_back(component_nodes.back());

        auto first_distances = dijkstra_from_sources(graph, {component_nodes.front()});
        size_t farthest = component_nodes.front();
        cumulative_time_int farthest_distance = 0;
        for (auto node : component_nodes) {
            if (first_distances[node] < INF_DISTANCE && first_distances[node] > farthest_distance) {
                farthest_distance = first_distances[node];
                farthest = node;
            }
        }
        landmark_candidates.push_back(farthest);

        cumulative_time_int best_upper_bound = INF_DISTANCE;
        std::sort(landmark_candidates.begin(), landmark_candidates.end());
        landmark_candidates.erase(std::unique(landmark_candidates.begin(), landmark_candidates.end()), landmark_candidates.end());
        for (auto landmark : landmark_candidates) {
            auto distances = landmark == component_nodes.front()
                                 ? first_distances
                                 : dijkstra_from_sources(graph, {static_cast<uint64_t>(landmark)});
            cumulative_time_int eccentricity = 0;
            bool finite = true;
            for (auto node : component_nodes) {
                auto distance = distances[node];
                if (distance >= INF_DISTANCE) {
                    finite = false;
                    break;
                }
                eccentricity = std::max(eccentricity, distance);
            }
            if (finite && eccentricity <= INF_DISTANCE / 2) {
                best_upper_bound = std::min(best_upper_bound, eccentricity * 2);
            }
        }
        graph_cache.interior_component_diameter_upper_bounds.push_back(best_upper_bound);
    }

    graph_cache.interior_diameter_upper_bound = 0;
    for (auto upper_bound : graph_cache.interior_component_diameter_upper_bounds) {
        if (upper_bound >= INF_DISTANCE) {
            graph_cache.interior_diameter_upper_bound = INF_DISTANCE;
            return;
        }
        graph_cache.interior_diameter_upper_bound = std::max(graph_cache.interior_diameter_upper_bound, upper_bound);
    }
}

std::optional<cumulative_time_int> component_graph_component_upper_bound(
    const ProcessingClusterGraphCache& graph_cache,
    const std::vector<uint64_t>& component_detectors) {
    if (component_detectors.empty() || graph_cache.interior_component_id_by_vertex.empty()) {
        return std::nullopt;
    }
    if (component_detectors[0] >= graph_cache.interior_component_id_by_vertex.size()) {
        return std::nullopt;
    }
    size_t component_id = graph_cache.interior_component_id_by_vertex[component_detectors[0]];
    if (component_id == SIZE_MAX || component_id >= graph_cache.interior_component_diameter_upper_bounds.size()) {
        return std::nullopt;
    }
    for (auto detector : component_detectors) {
        if (detector >= graph_cache.interior_component_id_by_vertex.size() ||
            graph_cache.interior_component_id_by_vertex[detector] != component_id) {
            return std::nullopt;
        }
    }
    return graph_cache.interior_component_diameter_upper_bounds[component_id];
}

bool clusters_influence_intersect(const ProcessingCluster& a_cluster, const ProcessingCluster& b_cluster) {
    if (a_cluster.influence_is_full_graph || b_cluster.influence_is_full_graph) {
        return true;
    }
    const auto& a = a_cluster.influence_vertices;
    const auto& b = b_cluster.influence_vertices;
    size_t ia = 0;
    size_t ib = 0;
    while (ia < a.size() && ib < b.size()) {
        if (a[ia] == b[ib]) {
            return true;
        }
        if (a[ia] < b[ib]) {
            ia++;
        } else {
            ib++;
        }
    }
    return false;
}

void assign_parent_pointers(ProcessingClusterHierarchy& hierarchy) {
    size_t full_root_id = NO_PROCESSING_CLUSTER_PARENT;
    size_t full_root_level = 0;
    for (const auto& cluster : hierarchy.clusters) {
        if (!cluster.influence_is_full_graph) {
            continue;
        }
        if (full_root_id != NO_PROCESSING_CLUSTER_PARENT) {
            full_root_id = NO_PROCESSING_CLUSTER_PARENT;
            break;
        }
        full_root_id = cluster.id;
        full_root_level = cluster.level;
    }
    if (full_root_id != NO_PROCESSING_CLUSTER_PARENT) {
        hierarchy.root_cluster_ids.clear();
        for (auto& cluster : hierarchy.clusters) {
            if (cluster.id != full_root_id && cluster.level < full_root_level) {
                cluster.parent_id = full_root_id;
            }
        }
        hierarchy.root_cluster_ids.push_back(full_root_id);
        return;
    }

    for (auto& cluster : hierarchy.clusters) {
        for (size_t level = cluster.level + 1; level < hierarchy.cluster_ids_by_level.size(); level++) {
            bool found = false;
            for (auto candidate_id : hierarchy.cluster_ids_by_level[level]) {
                const auto& candidate = hierarchy.clusters[candidate_id];
                if (clusters_influence_intersect(cluster, candidate)) {
                    cluster.parent_id = candidate_id;
                    found = true;
                    break;
                }
            }
            if (found) {
                break;
            }
        }
    }

    hierarchy.root_cluster_ids.clear();
    for (const auto& cluster : hierarchy.clusters) {
        if (cluster.parent_id == NO_PROCESSING_CLUSTER_PARENT) {
            hierarchy.root_cluster_ids.push_back(cluster.id);
        }
    }
}

void append_level_if_needed(ProcessingClusterHierarchy& hierarchy, size_t level) {
    while (hierarchy.cluster_ids_by_level.size() <= level) {
        hierarchy.cluster_ids_by_level.push_back({});
    }
}

void append_cluster(ProcessingClusterHierarchy& hierarchy, ProcessingCluster&& cluster) {
    append_level_if_needed(hierarchy, cluster.level);
    cluster.id = hierarchy.clusters.size();
    hierarchy.cluster_ids_by_level[cluster.level].push_back(cluster.id);
    hierarchy.clusters.push_back(std::move(cluster));
}


ProcessingClusterBounds processing_cluster_bounds_for_completion_level(
    const ProcessingClusterConfig& config,
    size_t effective_max_level) {
    size_t bounds_level = std::max<size_t>(1, effective_max_level);
    return processing_cluster_bounds_for_level(bounds_level, config);
}

void fill_forced_final_root_influence_region(
    const MatchingGraph& graph,
    const ProcessingClusterGraphCache* graph_cache,
    ProcessingCluster& cluster) {
    if (graph_cache != nullptr) {
        fill_full_graph_influence_region_from_cache(*graph_cache, cluster);
        return;
    }
    cluster.influence_is_full_graph = false;
    cluster.influence_vertices.resize(graph.nodes.size());
    std::iota(cluster.influence_vertices.begin(), cluster.influence_vertices.end(), 0);
    cluster.boundary_endpoint_vertices.clear();
    cluster.nearest_boundary_match_distance = nearest_boundary_match_distance_for_sources(
        graph, cluster.active_detectors, graph_cache);
}

void append_forced_single_final_root_cluster(
    const MatchingGraph& graph,
    ProcessingClusterHierarchy& hierarchy,
    const ProcessingClusterGraphCache* graph_cache,
    const ProcessingClusterConfig& config,
    size_t forced_root_level,
    std::vector<uint64_t>& residual,
    ProcessingClusterProfilingStats* profiling_stats) {
    auto bounds = processing_cluster_bounds_for_completion_level(
        config, effective_max_level_unvalidated(config));

    auto diameter_start = steady_clock::now();
    cumulative_time_int diameter = 0;
    if (graph_cache != nullptr && graph_cache->has_all_pairs_distances && graph_cache->all_interior_pairs_finite) {
        diameter = graph_cache->max_finite_interior_distance;
    } else if (graph_cache != nullptr && graph_cache->interior_diameter_upper_bound < INF_DISTANCE) {
        diameter = graph_cache->interior_diameter_upper_bound;
    } else {
        diameter = INF_DISTANCE;
    }
    if (profiling_stats != nullptr) {
        profiling_stats->diameter_check_wall_ns += elapsed_ns(diameter_start);
    }

    cumulative_time_int diameter_bound = bounds.diameter_bound;
    if (diameter < INF_DISTANCE) {
        diameter_bound = std::max(diameter_bound, diameter);
    }
    cumulative_time_int buffer_bound = bounds.buffer_bound;
    if (graph_cache != nullptr) {
        cumulative_time_int full_buffer_bound = INF_DISTANCE;
        if (graph_cache->interior_diameter_upper_bound < INF_DISTANCE &&
            graph_cache->max_boundary_edge_weight < INF_DISTANCE &&
            graph_cache->interior_diameter_upper_bound <= INF_DISTANCE - graph_cache->max_boundary_edge_weight) {
            full_buffer_bound = graph_cache->interior_diameter_upper_bound + graph_cache->max_boundary_edge_weight;
        } else if (graph_cache->has_all_pairs_distances && graph_cache->all_interior_pairs_finite &&
                   graph_cache->max_finite_interior_distance < INF_DISTANCE &&
                   graph_cache->max_boundary_edge_weight < INF_DISTANCE &&
                   graph_cache->max_finite_interior_distance <= INF_DISTANCE - graph_cache->max_boundary_edge_weight) {
            full_buffer_bound = graph_cache->max_finite_interior_distance + graph_cache->max_boundary_edge_weight;
        }
        if (full_buffer_bound < INF_DISTANCE) {
            buffer_bound = std::max(buffer_bound, full_buffer_bound);
        }
    }

    ProcessingCluster cluster{
        0,
        forced_root_level,
        std::move(residual),
        {},
        false,
        {},
        NO_PROCESSING_CLUSTER_PARENT,
        diameter_bound,
        buffer_bound,
        diameter,
        INF_DISTANCE,
        INF_DISTANCE,
        true,
    };
    auto influence_start = steady_clock::now();
    fill_forced_final_root_influence_region(graph, graph_cache, cluster);
    if (profiling_stats != nullptr) {
        profiling_stats->influence_region_wall_ns += elapsed_ns(influence_start);
    }
    append_cluster(hierarchy, std::move(cluster));
    residual.clear();
}

bool try_append_full_root_completion_cluster(
    ProcessingClusterHierarchy& hierarchy,
    const ProcessingClusterGraphCache* graph_cache,
    const ProcessingClusterConfig& config,
    size_t current_level,
    std::vector<uint64_t>& residual,
    ProcessingClusterProfilingStats* profiling_stats) {
    if (graph_cache == nullptr || residual.empty() || !graph_cache->has_radius_neighbor_precompute) {
        return false;
    }
    size_t target_level = graph_cache->radius_neighbor_max_level;
    size_t effective_max_level = effective_max_level_unvalidated(config);
    if (target_level <= current_level + 1 || target_level > effective_max_level) {
        return false;
    }
    auto bounds = processing_cluster_bounds_for_level(target_level, config);
    bool buffer_covers_root = can_use_full_graph_influence_fast_path(*graph_cache, bounds.buffer_bound) ||
        buffer_covers_all_interior_distances(*graph_cache, bounds.buffer_bound) ||
        buffer_covers_all_interior_distances_by_upper_bound(*graph_cache, bounds.buffer_bound);
    bool diameter_covers_root = diameter_bound_covers_all_interior_distances(*graph_cache, bounds.diameter_bound) ||
        (graph_cache->interior_diameter_upper_bound < INF_DISTANCE &&
         graph_cache->interior_diameter_upper_bound <= bounds.diameter_bound);
    if (!buffer_covers_root || !diameter_covers_root) {
        return false;
    }

    auto diameter_start = steady_clock::now();
    cumulative_time_int diameter = graph_cache->has_all_pairs_distances && graph_cache->all_interior_pairs_finite
        ? graph_cache->max_finite_interior_distance
        : graph_cache->interior_diameter_upper_bound;
    if (profiling_stats != nullptr) {
        profiling_stats->diameter_check_wall_ns += elapsed_ns(diameter_start);
    }
    if ((residual.size() & 1) != 0) {
        const auto boundary_inclusive_diameter =
            farthest_boundary_match_distance_for_sources_from_cache(*graph_cache, residual);
        if (boundary_inclusive_diameter > bounds.diameter_bound) {
            return false;
        }
    }

    ProcessingCluster cluster{
        0,
        target_level,
        std::move(residual),
        {},
        false,
        {},
        NO_PROCESSING_CLUSTER_PARENT,
        bounds.diameter_bound,
        bounds.buffer_bound,
        diameter,
        INF_DISTANCE,
        INF_DISTANCE,
        false,
    };
    auto influence_start = steady_clock::now();
    if (can_use_full_graph_influence_fast_path(*graph_cache, cluster.buffer_bound)) {
        fill_full_graph_influence_region_from_cache(*graph_cache, cluster);
    } else {
        fill_full_interior_influence_region_from_cache(*graph_cache, cluster);
    }
    if (profiling_stats != nullptr) {
        profiling_stats->influence_region_wall_ns += elapsed_ns(influence_start);
    }
    append_cluster(hierarchy, std::move(cluster));
    residual.clear();
    return true;
}

}  // namespace


ProcessingClusterBounds processing_cluster_bounds_for_level(size_t level, const ProcessingClusterConfig& config) {
    validate_config(config);
    return processing_cluster_bounds_for_level_unvalidated(level, config);
}

size_t processing_cluster_effective_max_level(const ProcessingClusterConfig& config) {
    validate_config(config);
    return effective_max_level_unvalidated(config);
}

bool processing_cluster_bounds_phi_buffer_certificate_ok(
    const ProcessingClusterBounds& bounds) {
    if (!bounds.has_phi_buffer_certificate || !(bounds.phi > 0) || !std::isfinite(bounds.phi)) {
        return false;
    }
    return bounds.buffer_bound >= bounds.phi_required_buffer_bound &&
           bounds.buffer_bound >= bounds.phi_budget_required_buffer_bound &&
           bounds.buffer_bound >= bounds.ratio_required_buffer_bound &&
           bounds.buffer_bound >= bounds.min_required_buffer_bound;
}

std::vector<double> processing_cluster_stopping_lemma_phi_values(
    const ProcessingClusterConfig& config) {
    if (effective_max_level_unvalidated(config) == 0) {
        return {};
    }
    size_t max_level = effective_max_level_unvalidated(config);
    std::vector<long double> phi_ld(max_level + 1, 0);
    std::vector<double> phi(max_level + 1, 0);
    phi_ld[1] = 1;
    phi[1] = 1;
    for (size_t level = 2; level <= max_level; level++) {
        long double occupied_fraction = 0;
        for (size_t lower = 1; lower < level; lower++) {
            auto bounds = processing_cluster_bounds_for_level_unvalidated(lower, config);
            long double d_plus_1 = certificate_diameter_bound(bounds) + 1;
            long double b_minus_1 = certificate_buffer_bound(bounds) - 1;
            long double denominator =
                (phi_ld[lower] + 1) * d_plus_1 + phi_ld[lower] * b_minus_1;
            if (!(denominator > 0)) {
                phi_ld[level] = -std::numeric_limits<long double>::infinity();
                phi[level] = -std::numeric_limits<double>::infinity();
                return phi;
            }
            occupied_fraction += ((phi_ld[lower] + 2) * d_plus_1) / denominator;
        }
        phi_ld[level] = 1 - occupied_fraction;
        phi[level] = static_cast<double>(phi_ld[level]);
    }
    return phi;
}

void validate_processing_cluster_stopping_lemma_bounds(const ProcessingClusterConfig& config) {
    size_t max_level = effective_max_level_unvalidated(config);
    if (max_level == 0) {
        throw std::invalid_argument(
            "ProcessingClusterConfig has no usable levels for stopping-lemma validation.");
    }

    std::vector<long double> phi(max_level + 1, 0);
    phi[1] = 1;
    ProcessingClusterBounds previous{0, 0};
    for (size_t level = 1; level <= max_level; level++) {
        auto bounds = processing_cluster_bounds_for_level_unvalidated(level, config);
        if (!(bounds.diameter_bound > 0)) {
            throw std::invalid_argument(
                "Stopping-lemma bounds require d_k > 0, but level " +
                std::to_string(level) + " has diameter_bound=" +
                std::to_string(bounds.diameter_bound) + ".");
        }
        if (!(bounds.buffer_bound > 1)) {
            throw std::invalid_argument(
                "Stopping-lemma bounds require b_k > 1, but level " +
                std::to_string(level) + " has buffer_bound=" +
                std::to_string(bounds.buffer_bound) + ".");
        }
        if (level > 1) {
            if (bounds.diameter_bound < previous.diameter_bound ||
                bounds.buffer_bound < previous.buffer_bound) {
                throw std::invalid_argument(
                    "Stopping-lemma bounds are expected to be monotone non-decreasing; level " +
                    std::to_string(level) + " is smaller than the previous level.");
            }
            long double occupied_fraction = 0;
            for (size_t lower = 1; lower < level; lower++) {
                auto lower_bounds = processing_cluster_bounds_for_level_unvalidated(lower, config);
                long double d_plus_1 = certificate_diameter_bound(lower_bounds) + 1;
                long double b_minus_1 = certificate_buffer_bound(lower_bounds) - 1;
                long double denominator =
                    (phi[lower] + 1) * d_plus_1 + phi[lower] * b_minus_1;
                if (!(denominator > 0)) {
                    throw std::invalid_argument(
                        "Stopping-lemma phi recursion has a non-positive denominator before level " +
                        std::to_string(level) + ".");
                }
                occupied_fraction += ((phi[lower] + 2) * d_plus_1) / denominator;
            }
            phi[level] = 1 - occupied_fraction;
        }
        long double ratio =
            2 * (certificate_diameter_bound(bounds) + 1) /
            (certificate_buffer_bound(bounds) - 1);
        if (!(phi[level] > ratio && ratio > 0)) {
            std::ostringstream oss;
            oss << "Stopping-lemma inequality failed at level " << level
                << ": phi_k=" << static_cast<double>(phi[level])
                << ", 2(d_k+1)/(b_k-1)=" << static_cast<double>(ratio)
                << ", d_k=" << bounds.diameter_bound
                << ", b_k=" << bounds.buffer_bound << ".";
            throw std::invalid_argument(oss.str());
        }
        previous = bounds;
    }
}

std::vector<ProcessingClusterBounds> make_processing_cluster_stopping_lemma_bounds(
    const std::vector<cumulative_time_int>& diameter_bounds,
    cumulative_time_int min_buffer_slack) {
    if (diameter_bounds.empty()) {
        return {};
    }
    if (min_buffer_slack < 1) {
        min_buffer_slack = 1;
    }
    std::vector<ProcessingClusterBounds> bounds;
    bounds.reserve(diameter_bounds.size());
    std::vector<long double> phi(diameter_bounds.size() + 1, 0);
    phi[1] = 1;
    cumulative_time_int previous_diameter = 0;
    cumulative_time_int previous_buffer = 0;
    for (size_t level = 1; level <= diameter_bounds.size(); level++) {
        cumulative_time_int d = diameter_bounds[level - 1];
        if (!(d > 0)) {
            throw std::invalid_argument(
                "make_processing_cluster_stopping_lemma_bounds requires every diameter bound to be positive.");
        }
        d = std::max(d, previous_diameter);
        if (level > 1) {
            long double occupied_fraction = 0;
            for (size_t lower = 1; lower < level; lower++) {
                auto lower_bounds = bounds[lower - 1];
                long double d_plus_1 = certificate_diameter_bound(lower_bounds) + 1;
                long double b_minus_1 = certificate_buffer_bound(lower_bounds) - 1;
                long double denominator =
                    (phi[lower] + 1) * d_plus_1 + phi[lower] * b_minus_1;
                if (!(denominator > 0)) {
                    throw std::invalid_argument(
                        "Cannot construct stopping-lemma bounds because phi recursion became invalid.");
                }
                occupied_fraction += ((phi[lower] + 2) * d_plus_1) / denominator;
            }
            phi[level] = 1 - occupied_fraction;
        }
        if (!(phi[level] > 0)) {
            throw std::invalid_argument(
                "Cannot construct stopping-lemma bounds because phi_k is non-positive at level " +
                std::to_string(level) + ".  Increase earlier buffer bounds or reduce the number of levels.");
        }
        long double strict_threshold =
            2.0L * (static_cast<long double>(d) + 1) / phi[level] + 1;
        auto max_distance = std::numeric_limits<cumulative_time_int>::max() / 8;
        if (strict_threshold >= static_cast<long double>(max_distance)) {
            throw std::invalid_argument(
                "Cannot construct stopping-lemma bounds without overflowing cumulative_time_int.");
        }
        cumulative_time_int phi_required_buffer =
            strict_threshold_to_distance_requirement(strict_threshold, min_buffer_slack);
        cumulative_time_int ratio_required_buffer = 0;
        cumulative_time_int min_required_buffer = d + 2;
        cumulative_time_int b = phi_required_buffer;
        b = std::max<cumulative_time_int>(b, min_required_buffer);
        b = std::max(b, previous_buffer);
        bounds.push_back(make_phi_certified_processing_cluster_bounds(
            d, b, phi[level], phi_required_buffer, ratio_required_buffer, min_required_buffer));
        previous_diameter = d;
        previous_buffer = b;
    }

    ProcessingClusterConfig check_config;
    check_config.max_level = bounds.size();
    check_config.explicit_bounds_by_level = bounds;
    validate_processing_cluster_stopping_lemma_bounds(check_config);
    return bounds;
}

std::vector<ProcessingClusterBounds> make_processing_cluster_linear_stopping_lemma_bounds(
    cumulative_time_int diameter_unit,
    size_t num_levels,
    cumulative_time_int min_buffer_slack) {
    if (!(diameter_unit > 0)) {
        throw std::invalid_argument(
            "make_processing_cluster_linear_stopping_lemma_bounds requires diameter_unit > 0.");
    }
    std::vector<cumulative_time_int> diameters;
    diameters.reserve(num_levels);
    auto max_distance = std::numeric_limits<cumulative_time_int>::max() / 8;
    for (size_t level = 1; level <= num_levels; level++) {
        if (level > static_cast<size_t>(max_distance / diameter_unit)) {
            throw std::invalid_argument(
                "make_processing_cluster_linear_stopping_lemma_bounds would overflow cumulative_time_int.");
        }
        diameters.push_back(static_cast<cumulative_time_int>(level * static_cast<size_t>(diameter_unit)));
    }
    return make_processing_cluster_stopping_lemma_bounds(diameters, min_buffer_slack);
}

std::vector<ProcessingClusterBounds> make_processing_cluster_disjoint_stopping_lemma_bounds_from_candidates(
    const std::vector<cumulative_time_int>& candidate_diameter_bounds,
    size_t max_levels,
    cumulative_time_int min_buffer_slack,
    double min_buffer_to_diameter_ratio,
    bool enforce_legacy_gap_filter,
    bool min_buffer_ratio_first_level_only,
    bool enforce_phi_floor_budget,
    double phi_floor,
    double phi_budget_per_level) {
    if (candidate_diameter_bounds.empty() || max_levels == 0) {
        return {};
    }
    if (min_buffer_slack < 1) {
        min_buffer_slack = 1;
    }
    if (!(min_buffer_to_diameter_ratio >= 1.0) || !std::isfinite(min_buffer_to_diameter_ratio)) {
        throw std::invalid_argument(
            "make_processing_cluster_disjoint_stopping_lemma_bounds_from_candidates requires "
            "min_buffer_to_diameter_ratio to be finite and at least 1.");
    }
    if (enforce_phi_floor_budget) {
        if (!(phi_floor > 0 && phi_floor < 1) || !std::isfinite(phi_floor)) {
            throw std::invalid_argument(
                "make_processing_cluster_disjoint_stopping_lemma_bounds_from_candidates requires "
                "phi_floor to be finite and in the open interval (0, 1) when phi-floor mode is enabled.");
        }
        if (!(phi_budget_per_level > 0) || !std::isfinite(phi_budget_per_level)) {
            throw std::invalid_argument(
                "make_processing_cluster_disjoint_stopping_lemma_bounds_from_candidates requires "
                "phi_budget_per_level to be finite and positive when phi-floor mode is enabled.");
        }
        if (!(phi_budget_per_level < 1.0 - phi_floor)) {
            throw std::invalid_argument(
                "make_processing_cluster_disjoint_stopping_lemma_bounds_from_candidates requires "
                "phi_budget_per_level < 1 - phi_floor so at least one future level can remain above the phi floor.");
        }
    }

    std::vector<cumulative_time_int> candidates;
    candidates.reserve(candidate_diameter_bounds.size());
    for (auto d : candidate_diameter_bounds) {
        if (d > 0) {
            candidates.push_back(d);
        }
    }
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    if (candidates.empty()) {
        return {};
    }

    auto max_distance = std::numeric_limits<cumulative_time_int>::max() / 8;
    std::vector<ProcessingClusterBounds> bounds;
    bounds.reserve(std::min(max_levels, candidates.size()));
    std::vector<long double> phi(max_levels + 1, 0);

    // Use the additive candidate column exactly as supplied, e.g.
    // d_1=L_origin, raw candidate d_{n+1}=d_n+L_origin.  If requested,
    // choose a subsequence from this raw additive column by applying the
    // legacy separation filter d_candidate >= 3 d_prev + 4 b_prev to the
    // previously selected bound.  The buffer bound is still chosen from the
    // stopping-lemma phi recursion.  The configured min_buffer_to_diameter_ratio
    // remains only a lower bound; the effective b_k/d_k ratio generally changes
    // with k through phi_k.  When min_buffer_ratio_first_level_only is enabled,
    // the ratio lower bound is applied only to b_1; for k>1, b_k is selected
    // from the phi recursion plus the basic d+2 and monotonicity constraints.
    // When enforce_phi_floor_budget is enabled, b_k is further raised so the
    // current level contributes at most phi_budget_per_level to future phi
    // consumption, and candidate selection stops before accepting a level with
    // phi_k < phi_floor.
    for (size_t candidate_index = 0;
         candidate_index < candidates.size() && bounds.size() < max_levels;
         candidate_index++) {
        size_t level = bounds.size() + 1;
        cumulative_time_int d = candidates[candidate_index];
        if (enforce_legacy_gap_filter && !bounds.empty()) {
            const auto& previous = bounds.back();
            long double required_next_d =
                3.0L * static_cast<long double>(previous.diameter_bound) +
                4.0L * static_cast<long double>(previous.buffer_bound);
            if (static_cast<long double>(d) < required_next_d) {
                continue;
            }
        }
        if (level == 1) {
            phi[level] = 1;
        } else {
            long double occupied_fraction = 0;
            for (size_t lower = 1; lower < level; lower++) {
                auto lower_bounds = bounds[lower - 1];
                long double d_plus_1 = certificate_diameter_bound(lower_bounds) + 1;
                long double b_minus_1 = certificate_buffer_bound(lower_bounds) - 1;
                long double denominator =
                    (phi[lower] + 1) * d_plus_1 + phi[lower] * b_minus_1;
                if (!(denominator > 0)) {
                    return bounds;
                }
                occupied_fraction += ((phi[lower] + 2) * d_plus_1) / denominator;
            }
            phi[level] = 1 - occupied_fraction;
        }
        if (!(phi[level] > 0)) {
            break;
        }
        if (enforce_phi_floor_budget && phi[level] < static_cast<long double>(phi_floor)) {
            break;
        }

        long double d_plus_1 = static_cast<long double>(d) + 1;
        long double strict_threshold = 2.0L * d_plus_1 / phi[level] + 1;
        cumulative_time_int phi_budget_required_buffer = 0;
        long double budget_threshold = 0;
        if (enforce_phi_floor_budget) {
            long double eta = static_cast<long double>(phi_budget_per_level);
            budget_threshold =
                1 + (d_plus_1 / phi[level]) *
                    (((phi[level] + 2) / eta) - (phi[level] + 1));
        }
        bool apply_ratio_lower_bound = !min_buffer_ratio_first_level_only || level == 1;
        long double ratio_threshold = apply_ratio_lower_bound
            ? std::ceil(min_buffer_to_diameter_ratio * static_cast<long double>(d))
            : 0;
        long double threshold = strict_threshold;
        if (apply_ratio_lower_bound) {
            threshold = std::max(threshold, ratio_threshold);
        }
        if (enforce_phi_floor_budget) {
            threshold = std::max(threshold, budget_threshold);
        }
        if (threshold >= static_cast<long double>(max_distance)) {
            break;
        }
        cumulative_time_int phi_required_buffer =
            strict_threshold_to_distance_requirement(strict_threshold, min_buffer_slack);
        if (enforce_phi_floor_budget) {
            if (budget_threshold >= static_cast<long double>(max_distance)) {
                break;
            }
            phi_budget_required_buffer =
                strict_threshold_to_distance_requirement(budget_threshold, min_buffer_slack);
        }
        cumulative_time_int ratio_required_buffer = apply_ratio_lower_bound
            ? ceil_long_double_to_distance(min_buffer_to_diameter_ratio * static_cast<long double>(d))
            : 0;
        cumulative_time_int min_required_buffer = d + 2;
        cumulative_time_int b = static_cast<cumulative_time_int>(std::floor(threshold)) + min_buffer_slack;
        b = std::max<cumulative_time_int>(b, phi_required_buffer);
        b = std::max<cumulative_time_int>(b, phi_budget_required_buffer);
        b = std::max<cumulative_time_int>(b, ratio_required_buffer);
        b = std::max<cumulative_time_int>(b, min_required_buffer);
        if (!bounds.empty()) {
            b = std::max(b, bounds.back().buffer_bound);
        }
        bounds.push_back(make_phi_certified_processing_cluster_bounds(
            d, b, phi[level], phi_required_buffer, ratio_required_buffer, min_required_buffer,
            phi_budget_required_buffer));
    }

    if (!bounds.empty()) {
        ProcessingClusterConfig check_config;
        check_config.max_level = bounds.size();
        check_config.explicit_bounds_by_level = bounds;
        validate_processing_cluster_stopping_lemma_bounds(check_config);
    }
    return bounds;
}


std::vector<ProcessingClusterBounds> make_processing_cluster_disjoint_stopping_lemma_bounds_direct_min_d(
    cumulative_time_int first_diameter_bound,
    size_t max_levels,
    cumulative_time_int min_buffer_slack,
    double min_buffer_to_diameter_ratio,
    bool min_buffer_ratio_first_level_only,
    bool enforce_phi_floor_budget,
    double phi_floor,
    double phi_budget_per_level) {
    if (max_levels == 0) {
        return {};
    }
    if (!(first_diameter_bound > 0)) {
        throw std::invalid_argument(
            "make_processing_cluster_disjoint_stopping_lemma_bounds_direct_min_d requires "
            "first_diameter_bound to be positive.");
    }
    if (min_buffer_slack < 1) {
        min_buffer_slack = 1;
    }
    if (!(min_buffer_to_diameter_ratio >= 1.0) || !std::isfinite(min_buffer_to_diameter_ratio)) {
        throw std::invalid_argument(
            "make_processing_cluster_disjoint_stopping_lemma_bounds_direct_min_d requires "
            "min_buffer_to_diameter_ratio to be finite and at least 1.");
    }
    if (enforce_phi_floor_budget) {
        if (!(phi_floor > 0 && phi_floor < 1) || !std::isfinite(phi_floor)) {
            throw std::invalid_argument(
                "make_processing_cluster_disjoint_stopping_lemma_bounds_direct_min_d requires "
                "phi_floor to be finite and in the open interval (0, 1) when phi-floor mode is enabled.");
        }
        if (!(phi_budget_per_level > 0) || !std::isfinite(phi_budget_per_level)) {
            throw std::invalid_argument(
                "make_processing_cluster_disjoint_stopping_lemma_bounds_direct_min_d requires "
                "phi_budget_per_level to be finite and positive when phi-floor mode is enabled.");
        }
        if (!(phi_budget_per_level < 1.0 - phi_floor)) {
            throw std::invalid_argument(
                "make_processing_cluster_disjoint_stopping_lemma_bounds_direct_min_d requires "
                "phi_budget_per_level < 1 - phi_floor so at least one future level can remain above the phi floor.");
        }
    }

    auto max_distance = std::numeric_limits<cumulative_time_int>::max() / 8;
    std::vector<ProcessingClusterBounds> bounds;
    bounds.reserve(max_levels);
    std::vector<long double> phi(max_levels + 1, 0);

    // Level 1 starts at first_diameter_bound.  Every later level is the exact
    // minimal integer allowed by the disjointness/gap precondition:
    //     d_{k+1} = 3 d_k + 4 b_k.
    // No additive candidate column is constructed or scanned.
    for (size_t level = 1; level <= max_levels; level++) {
        cumulative_time_int d = first_diameter_bound;
        if (level > 1) {
            const auto& previous = bounds.back();
            long double required_next_d =
                3.0L * static_cast<long double>(previous.diameter_bound) +
                4.0L * static_cast<long double>(previous.buffer_bound);
            if (!(required_next_d > 0) || required_next_d >= static_cast<long double>(max_distance)) {
                break;
            }
            d = ceil_long_double_to_distance(required_next_d);
        }

        if (level == 1) {
            phi[level] = 1;
        } else {
            long double occupied_fraction = 0;
            for (size_t lower = 1; lower < level; lower++) {
                auto lower_bounds = bounds[lower - 1];
                long double d_plus_1 = certificate_diameter_bound(lower_bounds) + 1;
                long double b_minus_1 = certificate_buffer_bound(lower_bounds) - 1;
                long double denominator =
                    (phi[lower] + 1) * d_plus_1 + phi[lower] * b_minus_1;
                if (!(denominator > 0)) {
                    return bounds;
                }
                occupied_fraction += ((phi[lower] + 2) * d_plus_1) / denominator;
            }
            phi[level] = 1 - occupied_fraction;
        }
        if (!(phi[level] > 0)) {
            break;
        }
        if (enforce_phi_floor_budget && phi[level] < static_cast<long double>(phi_floor)) {
            break;
        }

        long double d_plus_1 = static_cast<long double>(d) + 1;
        long double strict_threshold = 2.0L * d_plus_1 / phi[level] + 1;
        cumulative_time_int phi_budget_required_buffer = 0;
        long double budget_threshold = 0;
        if (enforce_phi_floor_budget) {
            long double eta = static_cast<long double>(phi_budget_per_level);
            budget_threshold =
                1 + (d_plus_1 / phi[level]) *
                    (((phi[level] + 2) / eta) - (phi[level] + 1));
        }
        bool apply_ratio_lower_bound = !min_buffer_ratio_first_level_only || level == 1;
        long double ratio_threshold = apply_ratio_lower_bound
            ? std::ceil(min_buffer_to_diameter_ratio * static_cast<long double>(d))
            : 0;
        long double threshold = strict_threshold;
        if (apply_ratio_lower_bound) {
            threshold = std::max(threshold, ratio_threshold);
        }
        if (enforce_phi_floor_budget) {
            threshold = std::max(threshold, budget_threshold);
        }
        if (threshold >= static_cast<long double>(max_distance)) {
            break;
        }
        cumulative_time_int phi_required_buffer =
            strict_threshold_to_distance_requirement(strict_threshold, min_buffer_slack);
        if (enforce_phi_floor_budget) {
            if (budget_threshold >= static_cast<long double>(max_distance)) {
                break;
            }
            phi_budget_required_buffer =
                strict_threshold_to_distance_requirement(budget_threshold, min_buffer_slack);
        }
        cumulative_time_int ratio_required_buffer = apply_ratio_lower_bound
            ? ceil_long_double_to_distance(min_buffer_to_diameter_ratio * static_cast<long double>(d))
            : 0;
        cumulative_time_int min_required_buffer = d + 2;
        cumulative_time_int b = static_cast<cumulative_time_int>(std::floor(threshold)) + min_buffer_slack;
        b = std::max<cumulative_time_int>(b, phi_required_buffer);
        b = std::max<cumulative_time_int>(b, phi_budget_required_buffer);
        b = std::max<cumulative_time_int>(b, ratio_required_buffer);
        b = std::max<cumulative_time_int>(b, min_required_buffer);
        if (!bounds.empty()) {
            b = std::max(b, bounds.back().buffer_bound);
        }
        bounds.push_back(make_phi_certified_processing_cluster_bounds(
            d, b, phi[level], phi_required_buffer, ratio_required_buffer, min_required_buffer,
            phi_budget_required_buffer));
    }

    if (!bounds.empty()) {
        ProcessingClusterConfig check_config;
        check_config.max_level = bounds.size();
        check_config.explicit_bounds_by_level = bounds;
        validate_processing_cluster_stopping_lemma_bounds(check_config);
    }
    return bounds;
}


std::vector<ProcessingClusterBounds> make_processing_cluster_disjoint_stopping_lemma_bounds_optimize_final_d(
    cumulative_time_int first_diameter_bound,
    size_t target_levels,
    cumulative_time_int min_buffer_slack,
    double min_buffer_to_diameter_ratio,
    bool min_buffer_ratio_first_level_only,
    double phi_floor,
    double phi_budget_total_fraction,
    size_t search_rounds) {
    if (target_levels == 0) {
        return {};
    }
    if (!(first_diameter_bound > 0)) {
        throw std::invalid_argument(
            "make_processing_cluster_disjoint_stopping_lemma_bounds_optimize_final_d requires "
            "first_diameter_bound to be positive.");
    }
    if (min_buffer_slack < 1) {
        min_buffer_slack = 1;
    }
    if (!(min_buffer_to_diameter_ratio >= 1.0) || !std::isfinite(min_buffer_to_diameter_ratio)) {
        throw std::invalid_argument(
            "make_processing_cluster_disjoint_stopping_lemma_bounds_optimize_final_d requires "
            "min_buffer_to_diameter_ratio to be finite and at least 1.");
    }
    if (!(phi_floor > 0 && phi_floor < 1) || !std::isfinite(phi_floor)) {
        throw std::invalid_argument(
            "make_processing_cluster_disjoint_stopping_lemma_bounds_optimize_final_d requires "
            "phi_floor to be finite and in the open interval (0, 1).");
    }
    if (!(phi_budget_total_fraction > 0 && phi_budget_total_fraction <= 1) ||
        !std::isfinite(phi_budget_total_fraction)) {
        throw std::invalid_argument(
            "make_processing_cluster_disjoint_stopping_lemma_bounds_optimize_final_d requires "
            "phi_budget_total_fraction to be finite and in the interval (0, 1].");
    }
    if (search_rounds == 0) {
        search_rounds = 1;
    }

    struct Evaluation {
        bool ok = false;
        cumulative_time_int final_diameter = 0;
        cumulative_time_int final_buffer = 0;
        long double sum_buffer = 0;
        std::vector<ProcessingClusterBounds> bounds;
    };

    auto max_distance = std::numeric_limits<cumulative_time_int>::max() / 8;
    const size_t budget_count = target_levels > 0 ? target_levels - 1 : 0;
    const long double total_budget =
        static_cast<long double>(1.0 - phi_floor) * static_cast<long double>(phi_budget_total_fraction);
    const long double eta_min = std::max(1e-18L, total_budget * 1e-10L /
        static_cast<long double>(std::max<size_t>(1, budget_count)));

    auto actual_phi_for_next_level = [](
        const std::vector<ProcessingClusterBounds>& bounds,
        const std::vector<long double>& phi,
        size_t level) -> long double {
        if (level == 1) {
            return 1;
        }
        long double occupied_fraction = 0;
        for (size_t lower = 1; lower < level; lower++) {
            const auto& lower_bounds = bounds[lower - 1];
            long double d_plus_1 = static_cast<long double>(lower_bounds.diameter_bound) + 1;
            long double b_minus_1 = static_cast<long double>(lower_bounds.buffer_bound) - 1;
            long double denominator =
                (phi[lower] + 1) * d_plus_1 + phi[lower] * b_minus_1;
            if (!(denominator > 0)) {
                return -std::numeric_limits<long double>::infinity();
            }
            occupied_fraction += ((phi[lower] + 2) * d_plus_1) / denominator;
        }
        return 1 - occupied_fraction;
    };

    auto evaluate = [&](const std::vector<long double>& eta) -> Evaluation {
        Evaluation result;
        if (eta.size() != budget_count) {
            return result;
        }
        for (long double e : eta) {
            if (!(e > 0) || !std::isfinite(static_cast<double>(e))) {
                return result;
            }
        }

        std::vector<ProcessingClusterBounds> bounds;
        bounds.reserve(target_levels);
        std::vector<long double> phi(target_levels + 1, 0);

        for (size_t level = 1; level <= target_levels; level++) {
            cumulative_time_int d = first_diameter_bound;
            if (level > 1) {
                const auto& previous = bounds.back();
                long double required_next_d =
                    3.0L * static_cast<long double>(previous.diameter_bound) +
                    4.0L * static_cast<long double>(previous.buffer_bound);
                if (!(required_next_d > 0) || required_next_d >= static_cast<long double>(max_distance)) {
                    return result;
                }
                d = ceil_long_double_to_distance(required_next_d);
            }

            phi[level] = actual_phi_for_next_level(bounds, phi, level);
            if (!(phi[level] > 0) || phi[level] + 1e-15L < static_cast<long double>(phi_floor)) {
                return result;
            }

            long double d_plus_1 = static_cast<long double>(d) + 1;
            long double strict_threshold = 2.0L * d_plus_1 / phi[level] + 1;
            cumulative_time_int phi_budget_required_buffer = 0;
            long double budget_threshold = 0;
            if (level < target_levels) {
                long double e = std::max(eta[level - 1], eta_min);
                budget_threshold =
                    1 + (d_plus_1 / phi[level]) *
                        (((phi[level] + 2) / e) - (phi[level] + 1));
            }
            bool apply_ratio_lower_bound = !min_buffer_ratio_first_level_only || level == 1;
            long double ratio_threshold = apply_ratio_lower_bound
                ? std::ceil(min_buffer_to_diameter_ratio * static_cast<long double>(d))
                : 0;
            long double threshold = strict_threshold;
            if (apply_ratio_lower_bound) {
                threshold = std::max(threshold, ratio_threshold);
            }
            if (level < target_levels) {
                threshold = std::max(threshold, budget_threshold);
            }
            if (!(threshold > 0) || threshold >= static_cast<long double>(max_distance)) {
                return result;
            }
            cumulative_time_int phi_required_buffer =
                strict_threshold_to_distance_requirement(strict_threshold, min_buffer_slack);
            if (level < target_levels) {
                if (budget_threshold >= static_cast<long double>(max_distance)) {
                    return result;
                }
                phi_budget_required_buffer =
                    strict_threshold_to_distance_requirement(budget_threshold, min_buffer_slack);
            }
            cumulative_time_int ratio_required_buffer = apply_ratio_lower_bound
                ? ceil_long_double_to_distance(min_buffer_to_diameter_ratio * static_cast<long double>(d))
                : 0;
            if (d > max_distance - 2) {
                return result;
            }
            cumulative_time_int min_required_buffer = d + 2;
            cumulative_time_int b = static_cast<cumulative_time_int>(std::floor(threshold)) + min_buffer_slack;
            b = std::max<cumulative_time_int>(b, phi_required_buffer);
            b = std::max<cumulative_time_int>(b, phi_budget_required_buffer);
            b = std::max<cumulative_time_int>(b, ratio_required_buffer);
            b = std::max<cumulative_time_int>(b, min_required_buffer);
            if (!bounds.empty()) {
                b = std::max(b, bounds.back().buffer_bound);
            }
            if (b >= max_distance) {
                return result;
            }
            result.sum_buffer += static_cast<long double>(b);
            bounds.push_back(make_phi_certified_processing_cluster_bounds(
                d, b, phi[level], phi_required_buffer, ratio_required_buffer, min_required_buffer,
                phi_budget_required_buffer));
        }

        ProcessingClusterConfig check_config;
        check_config.max_level = bounds.size();
        check_config.explicit_bounds_by_level = bounds;
        try {
            validate_processing_cluster_stopping_lemma_bounds(check_config);
        } catch (const std::exception&) {
            return result;
        }
        result.ok = bounds.size() == target_levels;
        if (result.ok) {
            result.final_diameter = bounds.back().diameter_bound;
            result.final_buffer = bounds.back().buffer_bound;
            result.bounds = std::move(bounds);
        }
        return result;
    };

    auto better = [](const Evaluation& a, const Evaluation& b) -> bool {
        if (!a.ok) {
            return false;
        }
        if (!b.ok) {
            return true;
        }
        if (a.final_diameter != b.final_diameter) {
            return a.final_diameter < b.final_diameter;
        }
        if (a.final_buffer != b.final_buffer) {
            return a.final_buffer < b.final_buffer;
        }
        return a.sum_buffer < b.sum_buffer;
    };

    if (budget_count == 0) {
        std::vector<long double> eta;
        Evaluation single = evaluate(eta);
        return single.ok ? single.bounds : std::vector<ProcessingClusterBounds>{};
    }

    const size_t variable_count = budget_count + 1;  // final entry is unused slack.
    auto make_state_from_weights = [&](const std::vector<long double>& weights, long double slack_fraction) {
        std::vector<long double> state(variable_count, 0);
        long double eta_total = total_budget * std::max(0.0L, std::min(1.0L, 1.0L - slack_fraction));
        long double weight_sum = 0;
        for (auto w : weights) {
            weight_sum += std::max(0.0L, w);
        }
        if (!(weight_sum > 0)) {
            weight_sum = static_cast<long double>(budget_count);
        }
        for (size_t i = 0; i < budget_count; i++) {
            long double w = i < weights.size() ? std::max(0.0L, weights[i]) : 1.0L;
            if (!(w > 0)) {
                w = 1.0L;
            }
            state[i] = std::max(eta_min, eta_total * w / weight_sum);
        }
        long double used = 0;
        for (size_t i = 0; i < budget_count; i++) {
            used += state[i];
        }
        if (used > total_budget) {
            long double scale = total_budget / used;
            used = 0;
            for (size_t i = 0; i < budget_count; i++) {
                state[i] *= scale;
                used += state[i];
            }
        }
        state[budget_count] = std::max(0.0L, total_budget - used);
        return state;
    };

    auto eta_from_state = [&](const std::vector<long double>& state) {
        std::vector<long double> eta(budget_count, eta_min);
        for (size_t i = 0; i < budget_count; i++) {
            eta[i] = std::max(eta_min, state[i]);
        }
        return eta;
    };

    auto evaluation_from_state = [&](const std::vector<long double>& state) {
        return evaluate(eta_from_state(state));
    };

    std::vector<std::vector<long double>> seeds;
    std::vector<long double> uniform(budget_count, 1.0L);
    seeds.push_back(make_state_from_weights(uniform, 0.0L));
    seeds.push_back(make_state_from_weights(uniform, 0.05L));

    std::vector<long double> increasing(budget_count, 1.0L);
    std::vector<long double> decreasing(budget_count, 1.0L);
    for (size_t i = 0; i < budget_count; i++) {
        increasing[i] = std::pow(2.0L, static_cast<long double>(i));
        decreasing[i] = std::pow(2.0L, static_cast<long double>(budget_count - 1 - i));
    }
    seeds.push_back(make_state_from_weights(increasing, 0.0L));
    seeds.push_back(make_state_from_weights(decreasing, 0.0L));
    seeds.push_back(make_state_from_weights(increasing, 0.05L));
    seeds.push_back(make_state_from_weights(decreasing, 0.05L));
    for (size_t i = 0; i < budget_count; i++) {
        std::vector<long double> focused(budget_count, 1e-3L);
        focused[i] = 1.0L;
        seeds.push_back(make_state_from_weights(focused, 0.0L));
        seeds.push_back(make_state_from_weights(focused, 0.10L));
    }

    const std::vector<long double> step_fractions = {
        0.75L, 0.50L, 0.35L, 0.25L, 0.15L, 0.10L,
        0.06L, 0.035L, 0.02L, 0.012L, 0.007L, 0.004L, 0.002L, 0.001L
    };

    Evaluation best;
    std::vector<long double> best_state;
    for (auto state : seeds) {
        Evaluation current = evaluation_from_state(state);
        if (better(current, best)) {
            best = current;
            best_state = state;
        }
        for (size_t round = 0; round < search_rounds; round++) {
            bool improved_this_round = false;
            for (long double step_fraction : step_fractions) {
                bool improved_this_step = true;
                size_t guard = 0;
                while (improved_this_step && guard < variable_count * variable_count * 2) {
                    guard++;
                    improved_this_step = false;
                    for (size_t from = 0; from < variable_count; from++) {
                        long double lower_bound = from < budget_count ? eta_min : 0.0L;
                        long double available = state[from] - lower_bound;
                        if (!(available > eta_min)) {
                            continue;
                        }
                        long double delta = available * step_fraction;
                        if (!(delta > eta_min)) {
                            continue;
                        }
                        for (size_t to = 0; to < variable_count; to++) {
                            if (from == to) {
                                continue;
                            }
                            std::vector<long double> candidate_state = state;
                            candidate_state[from] -= delta;
                            candidate_state[to] += delta;
                            Evaluation candidate = evaluation_from_state(candidate_state);
                            if (better(candidate, current)) {
                                state = std::move(candidate_state);
                                current = std::move(candidate);
                                improved_this_step = true;
                                improved_this_round = true;
                                if (better(current, best)) {
                                    best = current;
                                    best_state = state;
                                }
                                break;
                            }
                        }
                        if (improved_this_step) {
                            break;
                        }
                    }
                }
            }
            if (!improved_this_round) {
                break;
            }
        }
    }

    if (!best.ok) {
        return {};
    }

    // Re-evaluate the winning state once so the returned bounds reflect the
    // final rounded integer budget exactly, even if the best Evaluation moved.
    (void)best_state;
    return best.bounds;
}



std::vector<ProcessingClusterBounds> make_processing_cluster_disjoint_stopping_lemma_bounds_phi_sequence(
    cumulative_time_int first_diameter_bound,
    size_t max_levels,
    cumulative_time_int min_buffer_slack,
    double min_buffer_to_diameter_ratio,
    bool min_buffer_ratio_first_level_only,
    double phi2,
    double phi_floor,
    double phi_decay_q,
    long double diameter_growth_additive) {
    if (max_levels == 0) {
        return {};
    }
    if (!(first_diameter_bound > 0)) {
        throw std::invalid_argument(
            "make_processing_cluster_disjoint_stopping_lemma_bounds_phi_sequence requires "
            "first_diameter_bound to be positive.");
    }
    if (min_buffer_slack < 1) {
        min_buffer_slack = 1;
    }
    if (!(min_buffer_to_diameter_ratio >= 1.0) || !std::isfinite(min_buffer_to_diameter_ratio)) {
        throw std::invalid_argument(
            "make_processing_cluster_disjoint_stopping_lemma_bounds_phi_sequence requires "
            "min_buffer_to_diameter_ratio to be finite and at least 1.");
    }
    if (!(phi_floor > 0 && phi_floor < 1) || !std::isfinite(phi_floor)) {
        throw std::invalid_argument(
            "make_processing_cluster_disjoint_stopping_lemma_bounds_phi_sequence requires "
            "phi_floor to be finite and in the open interval (0, 1).");
    }
    if (!(phi2 > phi_floor && phi2 < 1) || !std::isfinite(phi2)) {
        throw std::invalid_argument(
            "make_processing_cluster_disjoint_stopping_lemma_bounds_phi_sequence requires "
            "phi2 to be finite and satisfy phi_floor < phi2 < 1.");
    }
    if (!(phi_decay_q > 0 && phi_decay_q < 1) || !std::isfinite(phi_decay_q)) {
        throw std::invalid_argument(
            "make_processing_cluster_disjoint_stopping_lemma_bounds_phi_sequence requires "
            "phi_decay_q to be finite and in the open interval (0, 1).");
    }

    std::vector<ProcessingClusterBounds> bounds;
    bounds.reserve(max_levels);
    std::vector<long double> phi(max_levels + 1, 0);
    std::vector<long double> schedule_d(max_levels + 1, 0);
    std::vector<long double> schedule_b(max_levels + 1, 0);

    auto target_phi_for_level = [&](size_t level) -> long double {
        if (level <= 1) {
            return 1.0L;
        }
        long double floor_ld = static_cast<long double>(phi_floor);
        long double phi2_ld = static_cast<long double>(phi2);
        long double q_ld = static_cast<long double>(phi_decay_q);
        return floor_ld + (phi2_ld - floor_ld) * std::pow(q_ld, static_cast<long double>(level - 2));
    };

    auto actual_phi_for_level = [&](size_t level) -> long double {
        if (level == 1) {
            return 1.0L;
        }
        long double occupied_fraction = 0;
        for (size_t lower = 1; lower < level; lower++) {
            long double d_plus_1 = schedule_d[lower] + 1;
            long double b_minus_1 = schedule_b[lower] - 1;
            long double denominator =
                (phi[lower] + 1) * d_plus_1 + phi[lower] * b_minus_1;
            if (!(denominator > 0) || !std::isfinite(denominator)) {
                return -std::numeric_limits<long double>::infinity();
            }
            long double term = ((phi[lower] + 2) * d_plus_1) / denominator;
            if (!std::isfinite(term)) {
                return -std::numeric_limits<long double>::infinity();
            }
            occupied_fraction += term;
        }
        return 1 - occupied_fraction;
    };

    auto rounded_requirement_ld = [&](long double threshold) -> long double {
        if (!(threshold > 0) || !std::isfinite(threshold)) {
            return std::numeric_limits<long double>::infinity();
        }
        return std::floor(threshold) + static_cast<long double>(min_buffer_slack);
    };

    for (size_t level = 1; level <= max_levels; level++) {
        long double d_ld = static_cast<long double>(first_diameter_bound);
        if (level > 1) {
            d_ld = 3.0L * schedule_d[level - 1] + 4.0L * schedule_b[level - 1] +
                diameter_growth_additive;
        }
        if (!(d_ld > 0) || !std::isfinite(d_ld)) {
            break;
        }
        schedule_d[level] = d_ld;

        phi[level] = actual_phi_for_level(level);
        long double target_phi = target_phi_for_level(level);
        if (!(phi[level] > 0) || !std::isfinite(phi[level]) || phi[level] + 1e-15L < target_phi) {
            break;
        }

        long double next_target_phi = target_phi_for_level(level + 1);
        long double eta = phi[level] - next_target_phi;
        if (!(eta > 0) || !std::isfinite(eta)) {
            break;
        }

        long double d_plus_1 = d_ld + 1;
        long double strict_threshold = 2.0L * d_plus_1 / phi[level] + 1;
        long double budget_threshold =
            1 + (d_plus_1 / phi[level]) *
                (((phi[level] + 2) / eta) - (phi[level] + 1));

        bool apply_ratio_lower_bound = !min_buffer_ratio_first_level_only || level == 1;
        long double ratio_threshold = apply_ratio_lower_bound
            ? std::ceil(static_cast<long double>(min_buffer_to_diameter_ratio) * d_ld)
            : 0;

        long double phi_required_buffer_ld = rounded_requirement_ld(strict_threshold);
        long double phi_budget_required_buffer_ld = rounded_requirement_ld(budget_threshold);
        long double ratio_required_buffer_ld = apply_ratio_lower_bound ? ratio_threshold : 0;
        long double min_required_buffer_ld = d_ld + 2;

        long double b_ld = phi_required_buffer_ld;
        b_ld = std::max(b_ld, phi_budget_required_buffer_ld);
        if (apply_ratio_lower_bound) {
            b_ld = std::max(b_ld, ratio_required_buffer_ld);
        }
        b_ld = std::max(b_ld, min_required_buffer_ld);
        if (level > 1) {
            b_ld = std::max(b_ld, schedule_b[level - 1]);
        }
        if (!(b_ld > 1) || !std::isfinite(b_ld)) {
            break;
        }
        schedule_b[level] = b_ld;

        cumulative_time_int d = ceil_long_double_to_distance(d_ld);
        cumulative_time_int b = ceil_long_double_to_distance(b_ld);
        cumulative_time_int phi_required_buffer = ceil_long_double_to_distance(phi_required_buffer_ld);
        cumulative_time_int phi_budget_required_buffer = ceil_long_double_to_distance(phi_budget_required_buffer_ld);
        cumulative_time_int ratio_required_buffer = apply_ratio_lower_bound
            ? ceil_long_double_to_distance(ratio_required_buffer_ld)
            : 0;
        cumulative_time_int min_required_buffer = ceil_long_double_to_distance(min_required_buffer_ld);

        bounds.push_back(make_phi_certified_processing_cluster_bounds(
            d, b, phi[level], phi_required_buffer, ratio_required_buffer, min_required_buffer,
            phi_budget_required_buffer, d_ld, b_ld));
    }

    if (!bounds.empty()) {
        ProcessingClusterConfig check_config;
        check_config.max_level = bounds.size();
        check_config.explicit_bounds_by_level = bounds;
        validate_processing_cluster_stopping_lemma_bounds(check_config);
    }
    return bounds;
}

std::vector<ProcessingClusterBounds> make_processing_cluster_parameter_schedule_bounds(
    cumulative_time_int first_diameter_bound,
    size_t max_levels,
    double phi_floor,
    double phi_decay_q) {
    if (!(phi_floor > 0 && phi_floor < 1) || !std::isfinite(phi_floor)) {
        throw std::invalid_argument(
            "make_processing_cluster_parameter_schedule_bounds requires phi_floor in (0, 1).");
    }
    if (!(phi_decay_q > 0 && phi_decay_q < 1) || !std::isfinite(phi_decay_q)) {
        throw std::invalid_argument(
            "make_processing_cluster_parameter_schedule_bounds requires phi_decay_q in (0, 1).");
    }
    const double phi2 = phi_floor + (1.0 - phi_floor) * phi_decay_q;
    auto bounds = make_processing_cluster_disjoint_stopping_lemma_bounds_phi_sequence(
        first_diameter_bound,
        max_levels,
        1,
        1.0,
        false,
        phi2,
        phi_floor,
        phi_decay_q,
        2.0L * (static_cast<long double>(first_diameter_bound) - 1.0L));
    for (size_t k = 0; k < bounds.size(); k++) {
        if (bounds[k].runtime_bounds_clamped) {
            bounds.resize(k + 1);
            break;
        }
    }
    return bounds;
}


std::vector<size_t> make_processing_cluster_even_detector_count_caps(size_t num_levels) {
    std::vector<size_t> caps;
    caps.reserve(num_levels);
    for (size_t level = 1; level <= num_levels; level++) {
        if (level > std::numeric_limits<size_t>::max() / 2) {
            throw std::invalid_argument(
                "make_processing_cluster_even_detector_count_caps would overflow size_t.");
        }
        caps.push_back(2 * level);
    }
    return caps;
}

size_t infer_processing_cluster_precompute_max_level(
    const ProcessingClusterGraphCache& graph_cache,
    const ProcessingClusterConfig& config) {
    validate_config(config);
    size_t configured_max = processing_cluster_effective_max_level(config);
    for (size_t level = 1; level <= configured_max; level++) {
        auto bounds = processing_cluster_bounds_for_level(level, config);
        bool buffer_covers_root = can_use_full_graph_influence_fast_path(graph_cache, bounds.buffer_bound) ||
            buffer_covers_all_interior_distances(graph_cache, bounds.buffer_bound) ||
            buffer_covers_all_interior_distances_by_upper_bound(graph_cache, bounds.buffer_bound);
        bool diameter_covers_root = diameter_bound_covers_all_interior_distances(graph_cache, bounds.diameter_bound) ||
            (graph_cache.interior_diameter_upper_bound < INF_DISTANCE &&
             graph_cache.interior_diameter_upper_bound <= bounds.diameter_bound);
        if (buffer_covers_root && diameter_covers_root) {
            return level;
        }
    }
    return configured_max;
}

ProcessingClusterRadiusNeighborList build_radius_neighbor_list_for_bound(
    const MatchingGraph& graph,
    const ProcessingClusterGraphCache& graph_cache,
    cumulative_time_int buffer_bound) {
    ProcessingClusterRadiusNeighborList result;
    result.buffer_bound = buffer_bound;
    result.offsets.assign(graph_cache.num_nodes + 1, 0);

    std::vector<std::vector<uint32_t>> rows(graph_cache.num_nodes);
    if (graph_cache.has_all_pairs_distances) {
        parallel_for_rows(graph_cache.num_nodes, [&](size_t source) {
            auto& row = rows[source];
            for (size_t destination = 0; destination < graph_cache.num_nodes; destination++) {
                if (destination == source) {
                    continue;
                }
                if (cached_distance(graph_cache, source, destination) <= buffer_bound) {
                    row.push_back(static_cast<uint32_t>(destination));
                }
            }
        });
    } else {
        parallel_for_rows(graph_cache.num_nodes, [&](size_t source) {
            auto distances = dijkstra_from_sources(graph, {static_cast<uint64_t>(source)}, buffer_bound);
            auto& row = rows[source];
            for (size_t destination = 0; destination < distances.size(); destination++) {
                if (destination == source) {
                    continue;
                }
                if (distances[destination] <= buffer_bound) {
                    row.push_back(static_cast<uint32_t>(destination));
                }
            }
        });
    }

    uint64_t total = 0;
    for (size_t source = 0; source < rows.size(); source++) {
        result.offsets[source] = total;
        total += static_cast<uint64_t>(rows[source].size());
    }
    result.offsets[rows.size()] = total;
    result.neighbors.reserve(static_cast<size_t>(total));
    for (auto& row : rows) {
        result.neighbors.insert(result.neighbors.end(), row.begin(), row.end());
    }
    return result;
}

void ensure_processing_cluster_radius_neighbors_precomputed(
    const MatchingGraph& graph,
    ProcessingClusterGraphCache& graph_cache,
    const ProcessingClusterConfig& config) {
    validate_config(config);
    if (const char* e = std::getenv("PYMATCHING_SKIP_RADIUS_NEIGHBOR_PRECOMPUTE")) {
        if (std::string(e) == "1") {
            graph_cache.radius_neighbor_max_level = 0;
            graph_cache.radius_neighbors_by_level.clear();
            graph_cache.radius_neighbor_bounds_by_level.clear();
            graph_cache.has_radius_neighbor_precompute = false;
            return;
        }
    }
    size_t max_level = infer_processing_cluster_precompute_max_level(graph_cache, config);
    if (graph_cache.has_radius_neighbor_precompute &&
        graph_cache.radius_neighbor_max_level >= max_level &&
        graph_cache.radius_neighbors_by_level.size() > max_level &&
        cached_bounds_match(graph_cache.radius_neighbor_bounds_by_level, config, max_level)) {
        return;
    }

    graph_cache.radius_neighbor_beta = config.beta;
    graph_cache.radius_neighbor_gamma = config.gamma;
    graph_cache.radius_neighbor_lambda = config.lambda;
    graph_cache.radius_neighbor_max_level = max_level;
    graph_cache.radius_neighbor_bounds_by_level =
        processing_cluster_bounds_vector_unvalidated(config, max_level);
    graph_cache.radius_neighbors_by_level.clear();
    graph_cache.radius_neighbors_by_level.resize(max_level + 1);

    if (graph_cache.has_all_pairs_distances) {
        graph_cache.radius_neighbors_by_level = build_radius_neighbor_lists_for_levels_from_all_pairs(
            graph_cache, config, max_level);
    } else {
        for (size_t level = 1; level <= max_level; level++) {
            auto bounds = processing_cluster_bounds_for_level(level, config);
            // If the buffer already covers the graph component by cheap graph-only
            // upper-bound metadata, the online clustering code can group by
            // interior_component_id directly.  Do not spend one-time precompute
            // building a huge near-complete radius-neighbor CSR for that level.
            if (can_use_full_graph_influence_fast_path(graph_cache, bounds.buffer_bound) ||
                buffer_covers_all_interior_distances_by_upper_bound(graph_cache, bounds.buffer_bound)) {
                graph_cache.radius_neighbors_by_level[level].buffer_bound = bounds.buffer_bound;
                graph_cache.radius_neighbors_by_level[level].offsets.assign(graph_cache.num_nodes + 1, 0);
                continue;
            }
            graph_cache.radius_neighbors_by_level[level] =
                build_radius_neighbor_list_for_bound(graph, graph_cache, bounds.buffer_bound);
        }
    }
    graph_cache.has_radius_neighbor_precompute = true;
}

ProcessingClusterGraphCache build_processing_cluster_graph_cache(const MatchingGraph& graph) {
    ProcessingClusterGraphCache graph_cache;
    graph_cache.num_nodes = graph.nodes.size();
    graph_cache.boundary_edge_weight_by_vertex.assign(graph.nodes.size(), INF_DISTANCE);
    graph_cache.nearest_boundary_match_distance_by_vertex.assign(graph.nodes.size(), INF_DISTANCE);

    uint64_t num_entries =
        (static_cast<uint64_t>(graph.nodes.size()) * static_cast<uint64_t>(graph.nodes.size() + 1)) / 2;
    const uint64_t cache_max_bytes = processing_cluster_cache_max_bytes();
    bool build_all_pairs = graph.nodes.size() <= processing_cluster_all_pairs_max_nodes() &&
        num_entries <= cache_max_bytes / sizeof(uint32_t);
    graph_cache.has_all_pairs_distances = build_all_pairs;
    if (build_all_pairs) {
        graph_cache.packed_all_pairs_interior_distances.assign(num_entries, PACKED_INF_DISTANCE);
        graph_cache.packed_all_pairs_row_offsets.resize(graph.nodes.size());
        for (size_t source = 0; source < graph.nodes.size(); source++) {
            graph_cache.packed_all_pairs_row_offsets[source] =
                static_cast<uint64_t>(source) * static_cast<uint64_t>(graph.nodes.size()) -
                (static_cast<uint64_t>(source) * static_cast<uint64_t>(source - 1)) / 2;
        }
    }

    for (size_t node = 0; node < graph.nodes.size(); node++) {
        for (size_t k = 0; k < graph.nodes[node].neighbors.size(); k++) {
            if (graph.nodes[node].neighbors[k] != nullptr) {
                continue;
            }
            auto edge_weight = static_cast<cumulative_time_int>(graph.nodes[node].neighbor_weights[k]);
            graph_cache.boundary_edge_weight_by_vertex[node] =
                std::min(graph_cache.boundary_edge_weight_by_vertex[node], edge_weight);
        }
        if (graph_cache.boundary_edge_weight_by_vertex[node] < INF_DISTANCE) {
            graph_cache.boundary_endpoint_vertices.push_back(node);
            graph_cache.max_boundary_edge_weight =
                std::max(graph_cache.max_boundary_edge_weight, graph_cache.boundary_edge_weight_by_vertex[node]);
        }
    }

    fill_interior_component_upper_bounds(graph, graph_cache);

    if (build_all_pairs) {
        for (size_t source = 0; source < graph.nodes.size(); source++) {
            auto distances = dijkstra_from_sources(graph, {source});
            for (size_t destination = 0; destination < graph.nodes.size(); destination++) {
                if (destination < source) {
                    continue;
                }
                auto distance = distances[destination];
                uint64_t index = graph_cache.packed_all_pairs_row_offsets[source] + (destination - source);
                if (distance >= INF_DISTANCE) {
                    graph_cache.all_interior_pairs_finite = false;
                    graph_cache.packed_all_pairs_interior_distances[index] = PACKED_INF_DISTANCE;
                } else if (distance > static_cast<cumulative_time_int>(PACKED_INF_DISTANCE - 1)) {
                    throw std::invalid_argument(
                        "Processing cluster graph cache cannot represent an interior distance larger than 32 bits.");
                } else {
                    graph_cache.max_finite_interior_distance =
                        std::max(graph_cache.max_finite_interior_distance, distance);
                    graph_cache.packed_all_pairs_interior_distances[index] = static_cast<uint32_t>(distance);
                }
            }
            for (auto endpoint : graph_cache.boundary_endpoint_vertices) {
                auto distance = distances[endpoint];
                auto boundary_weight = graph_cache.boundary_edge_weight_by_vertex[endpoint];
                if (distance >= INF_DISTANCE || boundary_weight >= INF_DISTANCE || distance > INF_DISTANCE - boundary_weight) {
                    continue;
                }
                graph_cache.nearest_boundary_match_distance_by_vertex[source] =
                    std::min(graph_cache.nearest_boundary_match_distance_by_vertex[source], distance + boundary_weight);
            }
        }
        graph_cache.interior_diameter_upper_bound = graph_cache.max_finite_interior_distance;
    } else {
        // Large graphs would spend minutes building an all-pairs cache. Keep the
        // cheap boundary-distance part of the cache so parity/boundary decisions
        // still avoid one Dijkstra per cluster, and let per-shot component
        // construction use the lightweight flooder.
        graph_cache.all_interior_pairs_finite = false;
        graph_cache.max_finite_interior_distance = INF_DISTANCE;
        graph_cache.nearest_boundary_match_distance_by_vertex = dijkstra_from_boundary_edges(graph);
    }
    return graph_cache;
}

ProcessingClusterHierarchy build_processing_cluster_hierarchy(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& active_detectors,
    const ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache,
    ProcessingClusterProfilingStats* profiling_stats) {
    validate_config(config);
    ProcessingClusterHierarchy hierarchy;
    size_t effective_max_level = effective_max_level_unvalidated(config);
    std::vector<uint64_t> residual = normalise_active_detectors(graph, active_detectors);
    std::optional<ActiveDetectorDistanceCache> active_distance_cache;
    auto ensure_active_distance_cache = [&]() -> ActiveDetectorDistanceCache* {
        if (active_distance_cache.has_value()) {
            return &*active_distance_cache;
        }
        if (graph_cache == nullptr || !graph_cache->has_all_pairs_distances || residual.empty() || residual.size() > 4096) {
            return nullptr;
        }
        active_distance_cache = build_active_detector_distance_cache(*graph_cache, residual);
        return &*active_distance_cache;
    };

    for (size_t level = 1; !residual.empty(); level++) {
        if (level > effective_max_level) {
            if (config.force_single_final_root) {
                std::sort(residual.begin(), residual.end());
                append_forced_single_final_root_cluster(
                    graph, hierarchy, graph_cache, config, effective_max_level + 1, residual, profiling_stats);
                break;
            }
            throw std::invalid_argument(
                "Failed to construct parity-aware processing clusters before the effective maximum level.");
        }
        auto bounds = processing_cluster_bounds_for_level(level, config);
        auto component_start = steady_clock::now();
        // The all-pairs cache wins for small/medium residual sets by avoiding a
        // multi-source Dijkstra, but it becomes quadratic in the number of active
        // detectors. Keep the graph flooder for larger residuals.
        FloodedThresholdComponents flooded_components;
        bool used_precomputable_radius_lookup_path = false;
        bool used_precomputed_radius_neighbor_list_model = false;
        uint64_t precomputable_radius_lookup_elapsed = 0;
        if (graph_cache != nullptr && (buffer_covers_all_interior_distances(*graph_cache, bounds.buffer_bound) || buffer_covers_all_interior_distances_by_upper_bound(*graph_cache, bounds.buffer_bound))) {
            // At high levels the buffer can cover the entire finite detector graph.
            // Then every remaining detection event is necessarily linked, so avoid
            // an O(active^2) cached-distance pass or a multi-source Dijkstra.
            flooded_components = build_single_threshold_component(residual.size());
        } else if (graph_cache != nullptr &&
                   radius_neighbors_for_level(*graph_cache, level, bounds.buffer_bound) != nullptr) {
            // The graph-only neighbor list was already materialized before the shot.
            // Prefer it in the canonical sparse-shot regime: level-1 rows are short,
            // while the all-pairs active scan is O(active^2).
            flooded_components = build_threshold_components_from_radius_neighbors(
                *graph_cache, *radius_neighbors_for_level(*graph_cache, level, bounds.buffer_bound), residual);
        } else if (graph_cache != nullptr && graph_cache->has_all_pairs_distances && residual.size() <= 512) {
            flooded_components = build_threshold_components_from_cache(*graph_cache, residual, bounds.buffer_bound);
        } else if (graph_cache != nullptr && graph_cache->has_all_pairs_distances && residual.size() <= 2048 && graph_cache->num_nodes <= 12000) {
            flooded_components = build_threshold_components_from_cache(*graph_cache, residual, bounds.buffer_bound);
        } else {
            // This branch is the online stand-in for querying a graph-only radius-neighbor
            // lookup at the current buffer threshold. The Dijkstra/flooding part depends
            // only on the detector graph metric and the threshold, and can be replaced by
            // code-distance-level preprocessing that materializes radius-neighbor lists.
            // Keep the paper-style level decomposition, but account this graph-metric
            // lookup construction separately so the ideal online path charges only the
            // active-set bookkeeping that remains after such preprocessing.
            auto radius_lookup_start = steady_clock::now();
            flooded_components = build_threshold_components_with_lightweight_cluster_flooder(graph, residual, bounds.buffer_bound);
            precomputable_radius_lookup_elapsed = elapsed_ns(radius_lookup_start);
            used_precomputable_radius_lookup_path = true;
        }
        if (profiling_stats != nullptr) {
            auto component_elapsed = elapsed_ns(component_start);
            profiling_stats->component_construction_wall_ns += component_elapsed;
            if (used_precomputable_radius_lookup_path) {
                profiling_stats->precomputable_distance_lookup_wall_ns += precomputable_radius_lookup_elapsed;
            } else if (used_precomputed_radius_neighbor_list_model) {
                // With a code-distance-level radius-neighbor lookup materialized for
                // this buffer threshold, component construction reduces to reading
                // the precomputed active-neighbor lists and unioning them. The current
                // benchmark still performs the equivalent cached-distance scan, so
                // charge that graph-metric thresholding work to preprocessing.
                profiling_stats->precomputable_distance_lookup_wall_ns += component_elapsed;
            }
        }
        std::vector<uint64_t> next_residual;
        next_residual.reserve(residual.size());

        for (const auto& component : flooded_components.components) {
            bool force_at_max_level = config.force_cluster_at_max_level && level >= effective_max_level;
            if (component.size() == 1 && !force_at_max_level) {
                uint64_t detector = residual[component[0]];
                std::vector<uint64_t> singleton_detectors{detector};
                const auto boundary_inclusive_diameter =
                    farthest_boundary_match_distance_for_sources(graph, singleton_detectors, graph_cache);
                if (boundary_inclusive_diameter > bounds.diameter_bound) {
                    next_residual.push_back(detector);
                    continue;
                }
                ProcessingCluster cluster{
                    0,
                    level,
                    std::move(singleton_detectors),
                    {},
                    false,
                    {},
                    NO_PROCESSING_CLUSTER_PARENT,
                    bounds.diameter_bound,
                    bounds.buffer_bound,
                    0,
                    INF_DISTANCE,
                    INF_DISTANCE,
                    false,
                };
                auto influence_start = steady_clock::now();
                if (graph_cache != nullptr &&
                    can_use_full_graph_influence_fast_path(*graph_cache, cluster.buffer_bound)) {
                    fill_full_graph_influence_region_from_cache(*graph_cache, cluster);
                } else if (graph_cache != nullptr &&
                    (buffer_covers_all_interior_distances(*graph_cache, cluster.buffer_bound) ||
                     buffer_covers_all_interior_distances_by_upper_bound(*graph_cache, cluster.buffer_bound))) {
                    fill_full_interior_influence_region_from_cache(*graph_cache, cluster);
                } else if (graph_cache != nullptr &&
                           radius_neighbors_for_level(*graph_cache, level, cluster.buffer_bound) != nullptr) {
                    fill_influence_region_from_radius_neighbors(
                        *graph_cache,
                        *radius_neighbors_for_level(*graph_cache, level, cluster.buffer_bound),
                        cluster);
                } else {
                    fill_influence_region(graph, cluster);
                }
                if (profiling_stats != nullptr) {
                    profiling_stats->influence_region_wall_ns += elapsed_ns(influence_start);
                }
                append_cluster(hierarchy, std::move(cluster));
                continue;
            }

            std::vector<uint64_t> component_detectors;
            component_detectors.reserve(component.size());
            for (auto index : component) {
                component_detectors.push_back(residual[index]);
            }
            size_t max_active_detectors = max_active_detectors_for_level_unvalidated(level, config);
            if (!force_at_max_level && max_active_detectors != 0 &&
                component_detectors.size() > max_active_detectors) {
                next_residual.insert(
                    next_residual.end(), component_detectors.begin(), component_detectors.end());
                continue;
            }
            auto diameter_start = steady_clock::now();
            cumulative_time_int diameter;
            bool used_online_precomputable_distance_lookup = false;
            if (graph_cache != nullptr && !force_at_max_level &&
                diameter_bound_covers_all_interior_distances(*graph_cache, bounds.diameter_bound)) {
                // Conservative upper bound is enough for the accept/reject test, and
                // avoids another O(component^2) distance pass for high-level clusters.
                diameter = graph_cache->max_finite_interior_distance;
            } else if (graph_cache != nullptr && !force_at_max_level &&
                graph_cache->interior_diameter_upper_bound <= bounds.diameter_bound) {
                // Landmark Dijkstras give a safe whole-graph diameter upper bound:
                // dist(u,v) <= dist(u,L)+dist(L,v). If the whole graph fits under
                // this level's diameter bound, then every component at this level
                // does too.
                diameter = graph_cache->interior_diameter_upper_bound;
            } else if (graph_cache != nullptr && !force_at_max_level) {
                auto component_upper_bound = component_graph_component_upper_bound(*graph_cache, component_detectors);
                if (component_upper_bound.has_value() && *component_upper_bound <= bounds.diameter_bound) {
                    diameter = *component_upper_bound;
                } else if (active_distance_cache.has_value()) {
                    diameter = component_diameter_from_active_distance_cache_with_cutoff(
                        *active_distance_cache,
                        component_detectors,
                        std::optional<cumulative_time_int>(bounds.diameter_bound));
                } else {
                    if (graph_cache == nullptr || !graph_cache->has_all_pairs_distances) {
                        used_online_precomputable_distance_lookup = true;
                        diameter = component_diameter_within_bound(graph, component_detectors, bounds.diameter_bound);
                    } else {
                        diameter = component_diameter_from_cache_with_cutoff(
                            *graph_cache,
                            component_detectors,
                            std::optional<cumulative_time_int>(bounds.diameter_bound));
                    }
                }
            } else if (active_distance_cache.has_value()) {
                diameter = component_diameter_from_active_distance_cache_with_cutoff(
                    *active_distance_cache,
                    component_detectors,
                    force_at_max_level ? std::nullopt
                                       : std::optional<cumulative_time_int>(bounds.diameter_bound));
            } else {
                if (graph_cache == nullptr || !graph_cache->has_all_pairs_distances) {
                    used_online_precomputable_distance_lookup = true;
                    diameter = component_diameter_within_bound(graph, component_detectors, bounds.diameter_bound);
                } else {
                    diameter = component_diameter_from_cache_with_cutoff(
                        *graph_cache,
                        component_detectors,
                        force_at_max_level ? std::nullopt
                                           : std::optional<cumulative_time_int>(bounds.diameter_bound));
                }
            }
            if (profiling_stats != nullptr) {
                auto diameter_elapsed = elapsed_ns(diameter_start);
                profiling_stats->diameter_check_wall_ns += diameter_elapsed;
                if (used_online_precomputable_distance_lookup) {
                    profiling_stats->precomputable_distance_lookup_wall_ns += diameter_elapsed;
                }
            }
            if (diameter <= bounds.diameter_bound || force_at_max_level) {
                if (!force_at_max_level && (component_detectors.size() & 1) != 0) {
                    const auto boundary_inclusive_diameter =
                        farthest_boundary_match_distance_for_sources(graph, component_detectors, graph_cache);
                    if (boundary_inclusive_diameter > bounds.diameter_bound) {
                        next_residual.insert(
                            next_residual.end(), component_detectors.begin(), component_detectors.end());
                        continue;
                    }
                }
                ProcessingCluster cluster{
                    0,
                    level,
                    std::move(component_detectors),
                    {},
                    false,
                    {},
                    NO_PROCESSING_CLUSTER_PARENT,
                    bounds.diameter_bound,
                    bounds.buffer_bound,
                    diameter,
                    INF_DISTANCE,
                    INF_DISTANCE,
                    force_at_max_level && diameter > bounds.diameter_bound,
                };
                auto influence_start = steady_clock::now();
                if (graph_cache != nullptr &&
                    can_use_full_graph_influence_fast_path(*graph_cache, cluster.buffer_bound)) {
                    fill_full_graph_influence_region_from_cache(*graph_cache, cluster);
                } else if (graph_cache != nullptr &&
                    (buffer_covers_all_interior_distances(*graph_cache, cluster.buffer_bound) ||
                     buffer_covers_all_interior_distances_by_upper_bound(*graph_cache, cluster.buffer_bound))) {
                    fill_full_interior_influence_region_from_cache(*graph_cache, cluster);
                } else if (graph_cache != nullptr &&
                           radius_neighbors_for_level(*graph_cache, level, cluster.buffer_bound) != nullptr) {
                    fill_influence_region_from_radius_neighbors(
                        *graph_cache,
                        *radius_neighbors_for_level(*graph_cache, level, cluster.buffer_bound),
                        cluster);
                } else {
                    fill_influence_region(graph, cluster);
                }
                if (profiling_stats != nullptr) {
                    profiling_stats->influence_region_wall_ns += elapsed_ns(influence_start);
                }
                append_cluster(hierarchy, std::move(cluster));
            } else {
                next_residual.insert(next_residual.end(), component_detectors.begin(), component_detectors.end());
            }
        }

        if (level >= effective_max_level && !next_residual.empty()) {
            if (config.force_single_final_root) {
                std::sort(next_residual.begin(), next_residual.end());
                append_forced_single_final_root_cluster(
                    graph, hierarchy, graph_cache, config, effective_max_level + 1, next_residual, profiling_stats);
                residual.clear();
                break;
            }
            throw std::invalid_argument(
                "Failed to construct parity-aware processing clusters before the effective maximum level.");
        }
        std::sort(next_residual.begin(), next_residual.end());
        if (try_append_full_root_completion_cluster(
                hierarchy, graph_cache, config, level, next_residual, profiling_stats)) {
            residual.clear();
            break;
        }
        residual = std::move(next_residual);
    }

    auto parent_assignment_start = steady_clock::now();
    assign_parent_pointers(hierarchy);
    if (config.force_single_final_root && hierarchy.root_cluster_ids.size() > 1) {
        std::vector<uint64_t> empty_final_root_sources;
        size_t forced_root_level = 1;
        for (const auto& cluster : hierarchy.clusters) {
            forced_root_level = std::max(forced_root_level, cluster.level + 1);
        }
        auto old_roots = hierarchy.root_cluster_ids;
        append_forced_single_final_root_cluster(
            graph, hierarchy, graph_cache, config, forced_root_level, empty_final_root_sources, profiling_stats);
        size_t final_root_id = hierarchy.clusters.empty() ? NO_PROCESSING_CLUSTER_PARENT : hierarchy.clusters.back().id;
        for (auto old_root_id : old_roots) {
            if (old_root_id < hierarchy.clusters.size() && old_root_id != final_root_id) {
                hierarchy.clusters[old_root_id].parent_id = final_root_id;
            }
        }
        hierarchy.root_cluster_ids.clear();
        if (final_root_id != NO_PROCESSING_CLUSTER_PARENT) {
            hierarchy.root_cluster_ids.push_back(final_root_id);
        }
    }
    if (profiling_stats != nullptr) {
        profiling_stats->parent_assignment_wall_ns += elapsed_ns(parent_assignment_start);
    }
    return hierarchy;
}

ProcessingClusterHierarchy build_root_direct_processing_cluster_hierarchy(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& active_detectors,
    const ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache,
    ProcessingClusterProfilingStats* profiling_stats) {
    validate_config(config);
    ProcessingClusterHierarchy candidate_hierarchy;
    size_t effective_max_level = effective_max_level_unvalidated(config);
    std::vector<uint64_t> all_active = normalise_active_detectors(graph, active_detectors);
    std::vector<uint64_t> residual = all_active;

    for (size_t level = config.root_direct_min_level; !residual.empty(); level++) {
        if (level > effective_max_level) {
            throw std::invalid_argument(
                "Failed to construct parity-aware root-direct processing clusters before the effective maximum level.");
        }
        auto bounds = processing_cluster_bounds_for_level(level, config);
        auto component_start = steady_clock::now();
        FloodedThresholdComponents flooded_components;
        bool used_precomputable_radius_lookup_path = false;
        bool used_precomputed_radius_neighbor_list_model = false;
        uint64_t precomputable_radius_lookup_elapsed = 0;
        if (graph_cache != nullptr && (buffer_covers_all_interior_distances(*graph_cache, bounds.buffer_bound) || buffer_covers_all_interior_distances_by_upper_bound(*graph_cache, bounds.buffer_bound))) {
            flooded_components = build_single_threshold_component(residual.size());
        } else if (graph_cache != nullptr &&
                   radius_neighbors_for_level(*graph_cache, level, bounds.buffer_bound) != nullptr) {
            flooded_components = build_threshold_components_from_radius_neighbors(
                *graph_cache, *radius_neighbors_for_level(*graph_cache, level, bounds.buffer_bound), residual);
        } else if (graph_cache != nullptr && graph_cache->has_all_pairs_distances && residual.size() <= 2048 && graph_cache->num_nodes <= 12000) {
            flooded_components = build_threshold_components_from_cache(*graph_cache, residual, bounds.buffer_bound);
            used_precomputed_radius_neighbor_list_model = true;
        } else {
            auto radius_lookup_start = steady_clock::now();
            flooded_components = build_threshold_components_with_lightweight_cluster_flooder(graph, residual, bounds.buffer_bound);
            precomputable_radius_lookup_elapsed = elapsed_ns(radius_lookup_start);
            used_precomputable_radius_lookup_path = true;
        }
        if (profiling_stats != nullptr) {
            auto component_elapsed = elapsed_ns(component_start);
            profiling_stats->component_construction_wall_ns += component_elapsed;
            if (used_precomputable_radius_lookup_path) {
                profiling_stats->precomputable_distance_lookup_wall_ns += precomputable_radius_lookup_elapsed;
            } else if (used_precomputed_radius_neighbor_list_model) {
                profiling_stats->precomputable_distance_lookup_wall_ns += component_elapsed;
            }
        }

        std::vector<uint64_t> next_residual;
        for (const auto& component : flooded_components.components) {
            std::vector<uint64_t> component_detectors;
            component_detectors.reserve(component.size());
            for (auto index : component) {
                component_detectors.push_back(residual[index]);
            }
            bool force_at_max_level = config.force_cluster_at_max_level && level >= effective_max_level;
            auto diameter_start = steady_clock::now();
            cumulative_time_int diameter = bounds.diameter_bound;
            bool used_online_precomputable_distance_lookup = false;
            if (!config.root_direct_skip_diameter_check) {
                if (graph_cache != nullptr && graph_cache->has_all_pairs_distances && !force_at_max_level &&
                    diameter_bound_covers_all_interior_distances(*graph_cache, bounds.diameter_bound)) {
                    diameter = graph_cache->max_finite_interior_distance;
                } else if (graph_cache != nullptr && !force_at_max_level &&
                    graph_cache->interior_diameter_upper_bound <= bounds.diameter_bound) {
                    diameter = graph_cache->interior_diameter_upper_bound;
                } else {
                    if (graph_cache == nullptr || !graph_cache->has_all_pairs_distances) {
                        used_online_precomputable_distance_lookup = true;
                        diameter = component_diameter_within_bound(graph, component_detectors, bounds.diameter_bound);
                    } else {
                        diameter = component_diameter_from_cache_with_cutoff(
                            *graph_cache,
                            component_detectors,
                            force_at_max_level ? std::nullopt
                                               : std::optional<cumulative_time_int>(bounds.diameter_bound));
                    }
                }
            }
            if (profiling_stats != nullptr) {
                auto diameter_elapsed = elapsed_ns(diameter_start);
                profiling_stats->diameter_check_wall_ns += diameter_elapsed;
                if (used_online_precomputable_distance_lookup) {
                    profiling_stats->precomputable_distance_lookup_wall_ns += diameter_elapsed;
                }
            }

            if (config.root_direct_skip_diameter_check || diameter <= bounds.diameter_bound || force_at_max_level) {
                cumulative_time_int nearest_boundary = INF_DISTANCE;
                cumulative_time_int boundary_inclusive_diameter = INF_DISTANCE;
                if (component_detectors.size() % 2 == 1) {
                    nearest_boundary = nearest_boundary_match_distance_for_sources(
                        graph, component_detectors, graph_cache);
                    boundary_inclusive_diameter = farthest_boundary_match_distance_for_sources(
                        graph, component_detectors, graph_cache);
                }
                if (cluster_can_stop_locally_from_parity_and_boundary(
                        component_detectors.size(), boundary_inclusive_diameter, bounds.diameter_bound)) {
                    ProcessingCluster cluster{
                        0,
                        level,
                        std::move(component_detectors),
                        {},
                        false,
                        {},
                        NO_PROCESSING_CLUSTER_PARENT,
                        bounds.diameter_bound,
                        bounds.buffer_bound,
                        diameter,
                        nearest_boundary,
                        INF_DISTANCE,
                        force_at_max_level && diameter > bounds.diameter_bound,
                    };
                    append_cluster(candidate_hierarchy, std::move(cluster));
                } else {
                    next_residual.insert(
                        next_residual.end(), component_detectors.begin(), component_detectors.end());
                }
            } else {
                next_residual.insert(next_residual.end(), component_detectors.begin(), component_detectors.end());
            }
        }

        if (level >= effective_max_level && !next_residual.empty()) {
            throw std::invalid_argument(
                "Failed to construct parity-aware root-direct processing clusters before the effective maximum level.");
        }
        std::sort(next_residual.begin(), next_residual.end());
        residual = std::move(next_residual);
    }

    ProcessingClusterHierarchy root_hierarchy;
    if (candidate_hierarchy.clusters.empty()) {
        return root_hierarchy;
    }

    std::vector<size_t> active_index_by_vertex(graph.nodes.size(), SIZE_MAX);
    for (size_t k = 0; k < all_active.size(); k++) {
        active_index_by_vertex[all_active[k]] = k;
    }
    std::vector<uint8_t> active_covered(all_active.size(), false);

    std::vector<size_t> candidate_ids(candidate_hierarchy.clusters.size());
    std::iota(candidate_ids.begin(), candidate_ids.end(), 0);
    std::stable_sort(candidate_ids.begin(), candidate_ids.end(), [&](size_t a, size_t b) {
        const auto& ca = candidate_hierarchy.clusters[a];
        const auto& cb = candidate_hierarchy.clusters[b];
        if (ca.level != cb.level) {
            return ca.level > cb.level;
        }
        return ca.active_detectors.size() > cb.active_detectors.size();
    });

    for (auto candidate_id : candidate_ids) {
        auto candidate = candidate_hierarchy.clusters[candidate_id];
        bool needed = false;
        for (auto detector : candidate.active_detectors) {
            if (detector < active_index_by_vertex.size() && !active_covered[active_index_by_vertex[detector]]) {
                needed = true;
                break;
            }
        }
        if (!needed) {
            continue;
        }

        auto influence_start = steady_clock::now();
        if (graph_cache != nullptr && can_use_full_graph_influence_fast_path(*graph_cache, candidate.buffer_bound)) {
            fill_full_graph_influence_region_from_cache(*graph_cache, candidate);
        } else {
            fill_influence_region(graph, candidate);
        }
        if (profiling_stats != nullptr) {
            profiling_stats->influence_region_wall_ns += elapsed_ns(influence_start);
        }
        if (candidate.influence_is_full_graph) {
            std::fill(active_covered.begin(), active_covered.end(), true);
        } else {
            for (auto vertex : candidate.influence_vertices) {
                if (vertex >= active_index_by_vertex.size()) {
                    continue;
                }
                auto active_index = active_index_by_vertex[vertex];
                if (active_index != SIZE_MAX) {
                    active_covered[active_index] = true;
                }
            }
        }
        append_cluster(root_hierarchy, std::move(candidate));
    }

    for (auto covered : active_covered) {
        if (!covered) {
            throw std::invalid_argument(
                "Root-direct processing cluster construction failed to cover every active detector.");
        }
    }
    root_hierarchy.root_cluster_ids.clear();
    for (const auto& cluster : root_hierarchy.clusters) {
        root_hierarchy.root_cluster_ids.push_back(cluster.id);
    }
    return root_hierarchy;
}

}  // namespace pm
