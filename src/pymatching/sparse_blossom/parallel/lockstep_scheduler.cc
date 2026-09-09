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

#include "pymatching/sparse_blossom/parallel/lockstep_scheduler.h"
#include "pymatching/sparse_blossom/parallel/growing_only_clusterer.h"

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <deque>
#include <functional>
#include <future>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <iostream>
#include <thread>
#include <string>
#include <unordered_map>
#if defined(__linux__)
#include <pthread.h>
#include <sched.h>
#endif

namespace pm {
namespace {

ClusterSubgraph build_worker_cluster_state(
    UserGraph& graph,
    const ProcessingCluster& cluster,
    const LockstepSchedulerConfig& config,
    bool allow_shared_global_mwpm = true);

using steady_clock = std::chrono::steady_clock;

template <typename T>
std::vector<T*> collect_live_lockstep_arena_objects(const Arena<T>& arena) {
    std::vector<T*> result = arena.allocated;
    std::sort(result.begin(), result.end());
    auto unused = arena.available;
    std::sort(unused.begin(), unused.end());
    result.erase(std::set_difference(
        result.begin(), result.end(),
        unused.begin(), unused.end(),
        result.begin()), result.end());
    return result;
}

size_t flooder_queue_bucket_event_count(const Mwpm& mwpm) {
    size_t total = 0;
    for (const auto& bucket : mwpm.flooder.queue.bit_buckets) {
        total += bucket.size();
    }
    return total;
}

size_t rebuild_flooder_queue_from_desired_trackers(Mwpm& mwpm) {
    auto& flooder = mwpm.flooder;
    flooder.queue.clear();
    size_t requeued = 0;
    for (auto& node : flooder.graph.nodes) {
        auto& tracker = node.node_event_tracker;
        if (!tracker.has_desired_time) {
            continue;
        }
        auto desired_time = tracker.desired_time;
        tracker.has_queued_time = false;
        tracker.set_desired_event(FloodCheckEvent(&node, desired_time), flooder.queue);
        requeued++;
    }
    auto live_regions = collect_live_lockstep_arena_objects(flooder.region_arena);
    for (auto* region : live_regions) {
        auto& tracker = region->shrink_event_tracker;
        if (!tracker.has_desired_time) {
            continue;
        }
        auto desired_time = tracker.desired_time;
        tracker.has_queued_time = false;
        tracker.set_desired_event(FloodCheckEvent(region, desired_time), flooder.queue);
        requeued++;
    }
    return requeued;
}

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

void maybe_pin_current_parentless_sparse_worker(size_t worker_index) {
#if defined(__linux__)
    const char* enabled = std::getenv("PYMATCHING_PARENTLESS_PIN_SPARSE_WORKERS");
    if (enabled == nullptr || std::string(enabled) != "1") {
        return;
    }
    size_t core_base = 0;
    if (const char* env = std::getenv("PYMATCHING_PARENTLESS_SPARSE_WORKER_CORE_BASE")) {
        core_base = static_cast<size_t>(std::stoull(env));
    }
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(static_cast<int>(core_base + worker_index), &cpuset);
    (void)pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
#else
    (void)worker_index;
#endif
}

class ParentlessLevelWorkerPool {
public:
    explicit ParentlessLevelWorkerPool(size_t worker_count) {
        if (worker_count < 2) {
            return;
        }
        workers.reserve(worker_count);
        for (size_t k = 0; k < worker_count; k++) {
            workers.emplace_back([this, k]() {
                maybe_pin_current_parentless_sparse_worker(k);
                worker_loop();
            });
        }
    }

    ~ParentlessLevelWorkerPool() {
        {
            std::lock_guard<std::mutex> lock(mutex);
            stopping = true;
        }
        cv.notify_all();
        for (auto& worker : workers) {
            if (worker.joinable()) {
                worker.join();
            }
        }
    }

    ParentlessLevelWorkerPool(const ParentlessLevelWorkerPool&) = delete;
    ParentlessLevelWorkerPool& operator=(const ParentlessLevelWorkerPool&) = delete;

    size_t size() const {
        return workers.size();
    }

    void parallel_for(size_t count, const std::function<void(size_t)>& body) {
        if (count == 0) {
            return;
        }
        if (workers.empty() || count == 1) {
            for (size_t k = 0; k < count; k++) {
                body(k);
            }
            return;
        }

        std::mutex done_mutex;
        std::condition_variable done_cv;
        std::exception_ptr first_exception = nullptr;
        size_t remaining = count;

        {
            std::lock_guard<std::mutex> lock(mutex);
            for (size_t k = 0; k < count; k++) {
                tasks.emplace_back([&, k]() {
                    try {
                        body(k);
                    } catch (...) {
                        std::lock_guard<std::mutex> done_lock(done_mutex);
                        if (!first_exception) {
                            first_exception = std::current_exception();
                        }
                    }
                    {
                        std::lock_guard<std::mutex> done_lock(done_mutex);
                        remaining--;
                        if (remaining == 0) {
                            done_cv.notify_one();
                        }
                    }
                });
            }
        }
        cv.notify_all();

        std::unique_lock<std::mutex> done_lock(done_mutex);
        done_cv.wait(done_lock, [&]() { return remaining == 0; });
        if (first_exception) {
            std::rethrow_exception(first_exception);
        }
    }

private:
    void worker_loop() {
        while (true) {
            std::function<void()> task;
            {
                std::unique_lock<std::mutex> lock(mutex);
                cv.wait(lock, [&]() { return stopping || !tasks.empty(); });
                if (stopping && tasks.empty()) {
                    return;
                }
                task = std::move(tasks.front());
                tasks.pop_front();
            }
            task();
        }
    }

    std::vector<std::thread> workers;
    std::deque<std::function<void()>> tasks;
    std::mutex mutex;
    std::condition_variable cv;
    bool stopping = false;
};

ParentlessLevelWorkerPool& parentless_level_worker_pool(size_t worker_count) {
    static std::mutex pool_mutex;
    static std::unique_ptr<ParentlessLevelWorkerPool> pool;
    static size_t pool_worker_count = 0;
    if (worker_count < 2) {
        static ParentlessLevelWorkerPool serial_pool(1);
        return serial_pool;
    }
    std::lock_guard<std::mutex> lock(pool_mutex);
    if (!pool || pool_worker_count != worker_count) {
        pool = std::make_unique<ParentlessLevelWorkerPool>(worker_count);
        pool_worker_count = worker_count;
    }
    return *pool;
}

ExtendedMatchingResult decode_full_graph(
    UserGraph& graph, const std::vector<uint64_t>& active_detectors, bool edge_correlations) {
    bool needs_search_graph = edge_correlations || graph.get_num_observables() > sizeof(pm::obs_int) * 8;
    auto& mwpm = needs_search_graph ? graph.get_mwpm_with_search_graph() : graph.get_mwpm();
    ExtendedMatchingResult result(graph.get_num_observables());
    decode_detection_events(mwpm, active_detectors, result.obs_crossed.data(), result.weight, edge_correlations);
    return result;
}

size_t find_maximum_level(const ProcessingClusterHierarchy& hierarchy) {
    if (hierarchy.cluster_ids_by_level.empty()) {
        return 0;
    }
    return hierarchy.cluster_ids_by_level.size() - 1;
}

size_t find_max_clusters_in_level(const ProcessingClusterHierarchy& hierarchy) {
    size_t best = 0;
    for (const auto& level : hierarchy.cluster_ids_by_level) {
        best = std::max(best, level.size());
    }
    return best;
}

size_t sum_cluster_active_detectors(const ProcessingClusterHierarchy& hierarchy) {
    size_t total = 0;
    for (const auto& cluster : hierarchy.clusters) {
        total += cluster.active_detectors.size();
    }
    return total;
}

size_t max_cluster_active_detectors(const ProcessingClusterHierarchy& hierarchy) {
    size_t best = 0;
    for (const auto& cluster : hierarchy.clusters) {
        best = std::max(best, cluster.active_detectors.size());
    }
    return best;
}


bool env_is_set_to_local(const char* name, const char* expected) {
    const char* value = std::getenv(name);
    return value != nullptr && std::string(value) == expected;
}

cumulative_time_int min_interior_edge_weight_for_direct_roots(const MatchingGraph& graph) {
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

cumulative_time_int cached_min_interior_edge_weight_for_direct_roots(const MatchingGraph& graph) {
    static thread_local std::unordered_map<const MatchingGraph*, cumulative_time_int> cache;
    auto it = cache.find(&graph);
    if (it != cache.end()) {
        return it->second;
    }
    auto value = min_interior_edge_weight_for_direct_roots(graph);
    cache.emplace(&graph, value);
    return value;
}

cumulative_time_int cached_min_nearest_boundary_for_direct_roots(
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

bool try_build_direct_forced_component_roots_hierarchy(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& active_detectors,
    const ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache,
    ProcessingClusterHierarchy& hierarchy,
    ProcessingClusterProfilingStats* profiling_stats) {
    if (config.adaptive_cluster_bounds ||
        env_is_set_to_local("PYMATCHING_GROWING_ONLY_DISABLE_DIRECT_FORCED_ROOTS", "1")) {
        return false;
    }
    if (env_is_set_to_local("PYMATCHING_GROWING_ONLY_KEEP_EMPTY_TERMINAL_ROOT", "1")) {
        return false;
    }
    if (env_is_set_to_local("PYMATCHING_GROWING_ONLY_FORCE_EXACT_AT_MAX_LEVEL", "1") ||
        env_is_set_to_local("PYMATCHING_SKIP_FORCED_MAX_EXACT_CLUSTER_METADATA", "0")) {
        return false;
    }
    if (!config.force_cluster_at_max_level || active_detectors.empty() || graph_cache == nullptr ||
        graph_cache->interior_component_id_by_vertex.size() != graph_cache->num_nodes) {
        return false;
    }
    const size_t max_level = processing_cluster_effective_max_level(config);
    if (max_level == 0) {
        return false;
    }

    const auto min_interior_edge = cached_min_interior_edge_weight_for_direct_roots(graph);
    const auto min_boundary = cached_min_nearest_boundary_for_direct_roots(graph_cache);
    for (size_t level = 1; level < max_level; level++) {
        const auto bounds = processing_cluster_bounds_for_level(level, config);
        if (!(bounds.buffer_bound < min_interior_edge && min_boundary > bounds.diameter_bound)) {
            return false;
        }
    }

    auto start = steady_clock::now();
    std::vector<uint64_t> residual = active_detectors;
    bool sorted_unique = true;
    for (size_t k = 0; k < residual.size(); k++) {
        auto detector = residual[k];
        if (detector >= graph.nodes.size()) {
            throw std::invalid_argument(
                "The detection event with index " + std::to_string(detector) +
                " does not correspond to a node in the graph, which only has " +
                std::to_string(graph.nodes.size()) + " nodes.");
        }
        if (detector < graph.is_user_graph_boundary_node.size() && graph.is_user_graph_boundary_node[detector]) {
            sorted_unique = false;
        }
        if (k != 0 && detector <= residual[k - 1]) {
            sorted_unique = false;
        }
    }
    if (!sorted_unique) {
        residual.erase(std::remove_if(
            residual.begin(), residual.end(), [&](uint64_t detector) {
                return detector < graph.is_user_graph_boundary_node.size() && graph.is_user_graph_boundary_node[detector];
            }), residual.end());
        std::sort(residual.begin(), residual.end());
        residual.erase(std::unique(residual.begin(), residual.end()), residual.end());
    }
    if (residual.empty()) {
        return false;
    }

    const auto root_bounds = processing_cluster_bounds_for_level(max_level, config);
    hierarchy.clusters.clear();
    hierarchy.cluster_ids_by_level.clear();
    hierarchy.root_cluster_ids.clear();
    hierarchy.cluster_ids_by_level.resize(max_level + 1);

    const size_t residual_size_for_stats = residual.size();
    std::vector<size_t> component_ids;
    std::vector<std::vector<uint64_t>> grouped_detectors;
    const bool forced_split_certificate_ok =
        processing_cluster_bounds_phi_buffer_certificate_ok(root_bounds);
    if (!forced_split_certificate_ok) {
        // Correctness fallback: when the phi/buffer certificate needed for
        // disjoint forced component roots is unavailable or fails, keep all
        // remaining active detectors in one trivial full-graph root.
        grouped_detectors.push_back(std::move(residual));
    } else {
        component_ids.reserve(4);
        grouped_detectors.reserve(4);
        for (auto detector_u64 : residual) {
            size_t detector = static_cast<size_t>(detector_u64);
            size_t component_id = detector < graph_cache->interior_component_id_by_vertex.size()
                                      ? graph_cache->interior_component_id_by_vertex[detector]
                                      : SIZE_MAX;
            size_t out = SIZE_MAX;
            if (component_id != SIZE_MAX) {
                for (size_t k = 0; k < component_ids.size(); k++) {
                    if (component_ids[k] == component_id) {
                        out = k;
                        break;
                    }
                }
            }
            if (out == SIZE_MAX) {
                component_ids.push_back(component_id);
                grouped_detectors.emplace_back();
                out = grouped_detectors.size() - 1;
            }
            grouped_detectors[out].push_back(detector_u64);
        }
    }

    size_t max_component_size = 0;
    for (auto& group : grouped_detectors) {
        if (group.empty()) {
            continue;
        }
        ProcessingCluster cluster;
        cluster.id = hierarchy.clusters.size();
        cluster.level = max_level;
        cluster.active_detectors = std::move(group);
        cluster.influence_is_full_graph = true;
        cluster.parent_id = NO_PROCESSING_CLUSTER_PARENT;
        cluster.diameter_bound = root_bounds.diameter_bound;
        cluster.buffer_bound = root_bounds.buffer_bound;
        cluster.diameter = 0;
        cluster.nearest_boundary_match_distance = std::numeric_limits<cumulative_time_int>::max();
        cluster.nearest_external_detector_distance = std::numeric_limits<cumulative_time_int>::max();
        cluster.forced_by_max_level = true;
        max_component_size = std::max(max_component_size, cluster.active_detectors.size());
        hierarchy.cluster_ids_by_level[max_level].push_back(cluster.id);
        hierarchy.root_cluster_ids.push_back(cluster.id);
        hierarchy.clusters.push_back(std::move(cluster));
    }
    if (hierarchy.clusters.empty()) {
        return false;
    }

    const auto elapsed = elapsed_ns(start);
    if (profiling_stats != nullptr) {
        // Keep accounting consistent with the existing component-lookup path:
        // candidate generation/grouping is reported as precomputable-lookup work
        // and component_construction_wall_ns remains the inclusive clustering cost.
        profiling_stats->component_construction_wall_ns += elapsed;
        profiling_stats->precomputable_distance_lookup_wall_ns += elapsed;
        profiling_stats->growing_only_used_component_lookup += 1;
        profiling_stats->growing_only_forced_max_level_shortcut_count += hierarchy.clusters.size();
        profiling_stats->growing_only_max_component_active_size = std::max<uint64_t>(
            profiling_stats->growing_only_max_component_active_size,
            static_cast<uint64_t>(max_component_size));
        profiling_stats->growing_only_max_residual_active_size = std::max<uint64_t>(
            profiling_stats->growing_only_max_residual_active_size,
            static_cast<uint64_t>(residual_size_for_stats));
    }
    return true;
}

ProcessingClusterHierarchy build_growing_only_processing_cluster_hierarchy(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& active_detectors,
    const ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache,
    ProcessingClusterProfilingStats* profiling_stats = nullptr) {
    ProcessingClusterHierarchy hierarchy;
    if (try_build_direct_forced_component_roots_hierarchy(
            graph, active_detectors, config, graph_cache, hierarchy, profiling_stats)) {
        return hierarchy;
    }

    std::vector<uint64_t> residual = active_detectors;
    std::sort(residual.begin(), residual.end());
    residual.erase(std::unique(residual.begin(), residual.end()), residual.end());

    ProcessingClusterConfig shot_cluster_config = config;
    const size_t max_level = processing_cluster_effective_max_level(shot_cluster_config);
    size_t highest_non_root_level = 0;

    auto append_cluster = [&](ProcessingCluster cluster) {
        cluster.id = hierarchy.clusters.size();
        if (hierarchy.cluster_ids_by_level.size() <= cluster.level) {
            hierarchy.cluster_ids_by_level.resize(cluster.level + 1);
        }
        hierarchy.cluster_ids_by_level[cluster.level].push_back(cluster.id);
        hierarchy.clusters.push_back(std::move(cluster));
        return hierarchy.clusters.back().id;
    };

    for (size_t level = 1; level <= max_level && !residual.empty(); level++) {
        auto level_start = steady_clock::now();
        auto result = run_growing_only_clustering_with_adaptive_bounds(
            graph, std::move(residual), level, shot_cluster_config, graph_cache);
        auto level_elapsed = elapsed_ns(level_start);
        if (profiling_stats != nullptr) {
            profiling_stats->component_construction_wall_ns += level_elapsed;
            profiling_stats->diameter_check_wall_ns += result.stats.event_processing_wall_ns;
            profiling_stats->precomputable_distance_lookup_wall_ns += result.stats.event_generation_wall_ns;
            profiling_stats->growing_only_used_active_pair_table += result.stats.used_active_pair_table ? 1 : 0;
            profiling_stats->growing_only_used_sparse_frontier += result.stats.used_sparse_frontier ? 1 : 0;
            profiling_stats->growing_only_used_component_lookup += result.stats.used_component_lookup ? 1 : 0;
            profiling_stats->growing_only_forced_max_level_shortcut_count += result.stats.forced_max_level_shortcut_count;
            profiling_stats->growing_only_processed_frontier_events += result.stats.processed_frontier_events;
            profiling_stats->growing_only_generated_collision_events += result.stats.generated_collision_events;
            profiling_stats->growing_only_processed_collision_events += result.stats.processed_collision_events;
            profiling_stats->growing_only_exact_internal_pair_checks += result.stats.exact_internal_pair_checks;
            profiling_stats->growing_only_exact_external_pair_checks += result.stats.exact_external_pair_checks;
            profiling_stats->growing_only_exact_diameter_checks += result.stats.exact_diameter_checks;
            profiling_stats->growing_only_max_component_active_size = std::max<uint64_t>(
                profiling_stats->growing_only_max_component_active_size,
                static_cast<uint64_t>(result.stats.max_component_active_size));
            profiling_stats->growing_only_max_residual_active_size = std::max<uint64_t>(
                profiling_stats->growing_only_max_residual_active_size,
                static_cast<uint64_t>(result.stats.max_residual_active_size));
        }
        const auto bounds = processing_cluster_bounds_for_level(level, config);
        for (auto& accepted : result.accepted_clusters) {
            ProcessingCluster cluster;
            cluster.level = level;
            cluster.active_detectors = std::move(accepted.active_detectors);
            cluster.influence_is_full_graph = true;
            cluster.parent_id = 0;  // patched to the final root id after it is appended
            cluster.diameter_bound = bounds.diameter_bound;
            cluster.buffer_bound = bounds.buffer_bound;
            cluster.diameter = accepted.diameter;
            cluster.nearest_boundary_match_distance = accepted.nearest_boundary_match_distance;
            cluster.nearest_external_detector_distance = accepted.first_external_collision_time >= GROWING_ONLY_INF_DISTANCE
                                                             ? std::numeric_limits<cumulative_time_int>::max()
                                                             : 2 * accepted.first_external_collision_time;
            cluster.forced_by_max_level = accepted.forced_by_max_level;
            append_cluster(std::move(cluster));
            highest_non_root_level = std::max(highest_non_root_level, level);
        }
        residual = std::move(result.residual_active_detectors);
    }

    if (hierarchy.clusters.empty()) {
        ProcessingCluster root;
        root.level = 0;
        root.active_detectors = residual;
        root.influence_is_full_graph = true;
        root.parent_id = NO_PROCESSING_CLUSTER_PARENT;
        root.diameter_bound = 0;
        root.buffer_bound = 0;
        root.diameter = 0;
        append_cluster(std::move(root));
        hierarchy.root_cluster_ids.push_back(0);
        return hierarchy;
    }

    // If growing-only clustering consumed the entire syndrome, do not append an
    // empty terminal root.  The previous hierarchy added an active-detector-free
    // parent and then imported all finished max-level children into it, which was
    // pure bookkeeping in the full-graph parentless pipeline.  Promote the final
    // accepted components to roots instead.  The root aggregation path already
    // supports multiple roots by summing compact/extended results, and this keeps
    // the old behavior available for debugging.
    const char* keep_empty_root_env = std::getenv("PYMATCHING_GROWING_ONLY_KEEP_EMPTY_TERMINAL_ROOT");
    const bool keep_empty_root = keep_empty_root_env != nullptr && std::string(keep_empty_root_env) == "1";
    if (residual.empty() && !keep_empty_root) {
        for (auto& cluster : hierarchy.clusters) {
            cluster.parent_id = NO_PROCESSING_CLUSTER_PARENT;
            hierarchy.root_cluster_ids.push_back(cluster.id);
        }
        return hierarchy;
    }

    ProcessingCluster root;
    root.level = std::max<size_t>(highest_non_root_level + 1, 1);
    root.active_detectors = residual;
    root.influence_is_full_graph = true;
    root.parent_id = NO_PROCESSING_CLUSTER_PARENT;
    auto root_bounds = processing_cluster_bounds_for_level(
        std::min(root.level, std::max<size_t>(max_level, 1)), shot_cluster_config);
    root.diameter_bound = root_bounds.diameter_bound;
    root.buffer_bound = root_bounds.buffer_bound;
    root.diameter = 0;
    const size_t root_id = append_cluster(std::move(root));
    for (size_t k = 0; k < hierarchy.clusters.size(); k++) {
        if (k == root_id) {
            continue;
        }
        hierarchy.clusters[k].parent_id = root_id;
    }
    hierarchy.root_cluster_ids.push_back(root_id);
    return hierarchy;
}

void finalize_cluster_step_parallel_accounting(
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config) {
    auto& stats = result.profiling_stats;
    stats.cluster_step_wall_ns = 0;
    stats.max_cluster_step_wall_ns = 0;
    stats.max_level_cluster_step_wall_ns = 0;
    stats.max_level_cluster_lifecycle_wall_ns = 0;
    stats.max_level_cluster_initial_setup_wall_ns = 0;
    stats.max_level_cluster_child_import_wall_ns = 0;
    stats.max_level_cluster_checkpoint_capture_wall_ns = 0;
    stats.max_level_cluster_checkpoint_restore_wall_ns = 0;
    stats.max_level_cluster_mark_stop_wall_ns = 0;
    stats.max_level_cluster_id = 0;
    stats.ideal_worker_cluster_makespan_wall_ns = 0;
    stats.parallel_cluster_events = 0;
    stats.parallel_level_max_event_sum = 0;
    stats.parallel_max_cluster_events = 0;
    stats.parallel_nonempty_event_levels = 0;
    stats.parallel_level_critical_path_algorithmic_time = 0;

    const auto cluster_count = result.hierarchy.clusters.size();
    if (result.cluster_step_wall_ns.size() < cluster_count) result.cluster_step_wall_ns.resize(cluster_count, 0);
    if (result.cluster_step_event_counts.size() < cluster_count) result.cluster_step_event_counts.resize(cluster_count, 0);
    if (result.cluster_critical_path_event_counts.size() < cluster_count) result.cluster_critical_path_event_counts.resize(cluster_count, 0);
    if (result.cluster_has_critical_path_event_count.size() < cluster_count) result.cluster_has_critical_path_event_count.resize(cluster_count, 0);
    if (result.cluster_has_cutoff_algorithmic_time.size() < cluster_count) result.cluster_has_cutoff_algorithmic_time.resize(cluster_count, 0);
    if (result.cluster_cutoff_algorithmic_times.size() < cluster_count) result.cluster_cutoff_algorithmic_times.resize(cluster_count, 0);
    if (result.cluster_prior_level_critical_event_counts.size() < cluster_count) result.cluster_prior_level_critical_event_counts.resize(cluster_count, 0);
    if (result.cluster_overlap_event_bases.size() < cluster_count) result.cluster_overlap_event_bases.resize(cluster_count, 0);
    if (result.cluster_events_at_or_before_cutoff.size() < cluster_count) result.cluster_events_at_or_before_cutoff.resize(cluster_count, 0);
    if (result.cluster_events_after_cutoff.size() < cluster_count) result.cluster_events_after_cutoff.resize(cluster_count, 0);
    if (result.cluster_stop_algorithmic_times.size() < cluster_count) result.cluster_stop_algorithmic_times.resize(cluster_count, 0);
    if (result.cluster_critical_path_algorithmic_times.size() < cluster_count) result.cluster_critical_path_algorithmic_times.resize(cluster_count, 0);
    if (result.cluster_initial_setup_wall_ns.size() < cluster_count) result.cluster_initial_setup_wall_ns.resize(cluster_count, 0);
    if (result.cluster_child_import_wall_ns.size() < cluster_count) result.cluster_child_import_wall_ns.resize(cluster_count, 0);
    if (result.cluster_checkpoint_capture_wall_ns.size() < cluster_count) result.cluster_checkpoint_capture_wall_ns.resize(cluster_count, 0);
    if (result.cluster_checkpoint_restore_wall_ns.size() < cluster_count) result.cluster_checkpoint_restore_wall_ns.resize(cluster_count, 0);
    if (result.cluster_mark_stop_wall_ns.size() < cluster_count) result.cluster_mark_stop_wall_ns.resize(cluster_count, 0);
    size_t max_level = 0;
    for (const auto& cluster : result.hierarchy.clusters) {
        max_level = std::max(max_level, cluster.level);
    }
    std::vector<uint64_t> all_cluster_work;
    all_cluster_work.reserve(cluster_count);
    std::vector<uint8_t> level_has_events(max_level + 1, 0);
    for (size_t cluster_id = 0; cluster_id < cluster_count; cluster_id++) {
        const uint64_t work_ns = result.cluster_step_wall_ns[cluster_id];
        const size_t level = result.hierarchy.clusters[cluster_id].level;
        const uint64_t raw_event_count = result.cluster_step_event_counts[cluster_id];
        const uint64_t critical_event_count = result.cluster_critical_path_event_counts[cluster_id];
        const bool has_critical_count = result.cluster_has_critical_path_event_count[cluster_id] != 0;
        const uint64_t event_count_for_max = has_critical_count ? critical_event_count : raw_event_count;
        stats.parallel_cluster_events += raw_event_count;
        stats.parallel_max_cluster_events = std::max(stats.parallel_max_cluster_events, event_count_for_max);
        stats.parallel_level_max_event_sum = std::max(stats.parallel_level_max_event_sum, event_count_for_max);
        stats.parallel_level_critical_path_algorithmic_time = std::max(
            stats.parallel_level_critical_path_algorithmic_time,
            result.cluster_critical_path_algorithmic_times[cluster_id]);
        if (raw_event_count > 0) {
            level_has_events[level] = 1;
        }
        stats.cluster_step_wall_ns += work_ns;
        stats.max_cluster_step_wall_ns = std::max(stats.max_cluster_step_wall_ns, work_ns);
        // Historical field name, now deliberately used as the all-level critical
        // path cluster step time. The measurement policy is all-level/all-cluster
        // overlap, so we take max over every cluster, not only max-level clusters.
        stats.max_level_cluster_step_wall_ns = std::max(stats.max_level_cluster_step_wall_ns, work_ns);

        // Include child imports and checkpoint/export work because they happen
        // during sparse-blossom execution, but exclude initial MWPM setup and
        // graph distribution. Take the slowest cluster over all levels.
        const uint64_t lifecycle_ns =
            result.cluster_step_wall_ns[cluster_id] +
            result.cluster_child_import_wall_ns[cluster_id] +
            result.cluster_checkpoint_capture_wall_ns[cluster_id] +
            result.cluster_checkpoint_restore_wall_ns[cluster_id] +
            result.cluster_mark_stop_wall_ns[cluster_id];
        if (lifecycle_ns >= stats.max_level_cluster_lifecycle_wall_ns) {
            stats.max_level_cluster_lifecycle_wall_ns = lifecycle_ns;
            stats.max_level_cluster_initial_setup_wall_ns = result.cluster_initial_setup_wall_ns[cluster_id];
            stats.max_level_cluster_child_import_wall_ns = result.cluster_child_import_wall_ns[cluster_id];
            stats.max_level_cluster_checkpoint_capture_wall_ns = result.cluster_checkpoint_capture_wall_ns[cluster_id];
            stats.max_level_cluster_checkpoint_restore_wall_ns = result.cluster_checkpoint_restore_wall_ns[cluster_id];
            stats.max_level_cluster_mark_stop_wall_ns = result.cluster_mark_stop_wall_ns[cluster_id];
            stats.max_level_cluster_id = cluster_id;
        }
        all_cluster_work.push_back(work_ns);
    }

    // parallel_level_max_event_sum is already the maximum reported event count.
    // For the parentless level pipeline, this is the final level-by-level
    // critical-path count. For other paths it falls back to the all-cluster max.
    for (size_t level = 0; level < level_has_events.size(); level++) {
        stats.parallel_nonempty_event_levels += level_has_events[level] ? 1 : 0;
    }

    size_t workers = config.ideal_cluster_worker_count ? config.ideal_cluster_worker_count : config.num_workers;
    workers = std::max<size_t>(1, workers);
    stats.ideal_cluster_worker_count_used = workers;
    if (!all_cluster_work.empty()) {
        std::sort(all_cluster_work.begin(), all_cluster_work.end(), std::greater<uint64_t>());
        const size_t used_workers = std::min(workers, all_cluster_work.size());
        std::vector<uint64_t> loads(used_workers, 0);
        for (uint64_t work_ns : all_cluster_work) {
            auto it = std::min_element(loads.begin(), loads.end());
            *it += work_ns;
        }
        stats.ideal_worker_cluster_makespan_wall_ns = *std::max_element(loads.begin(), loads.end());
    }
}

void expand_parent_influence_to_cover_direct_children(ProcessingClusterHierarchy& hierarchy, size_t num_graph_nodes) {
    // The stopped-state import code remaps every detector node owned by a child
    // stopped state into the parent subgraph.  The fast clustering code can keep
    // parent influence regions minimal; before building subgraphs, explicitly
    // make each parent subgraph contain the influence support of its children.
    // Children are constructed before parents in the current hierarchy builder.
    // Iterate in increasing cluster id order so that child influence is first
    // expanded by its descendants before the child itself is propagated into its
    // parent. This makes parent support transitive, which is required when a
    // stopped child state already contains imported descendant state.
    for (size_t child_id = 0; child_id < hierarchy.clusters.size(); child_id++) {
        auto parent_id = hierarchy.clusters[child_id].parent_id;
        if (parent_id == NO_PROCESSING_CLUSTER_PARENT || parent_id >= hierarchy.clusters.size()) {
            continue;
        }
        auto& parent = hierarchy.clusters[parent_id];
        const auto& child = hierarchy.clusters[child_id];
        if (parent.influence_is_full_graph) {
            continue;
        }
        if (child.influence_is_full_graph) {
            parent.influence_is_full_graph = true;
            parent.influence_vertices.clear();
            continue;
        }
        parent.influence_vertices.insert(
            parent.influence_vertices.end(), child.influence_vertices.begin(), child.influence_vertices.end());
        parent.influence_vertices.insert(
            parent.influence_vertices.end(), child.boundary_endpoint_vertices.begin(), child.boundary_endpoint_vertices.end());
        parent.influence_vertices.insert(
            parent.influence_vertices.end(), child.active_detectors.begin(), child.active_detectors.end());
        std::sort(parent.influence_vertices.begin(), parent.influence_vertices.end());
        parent.influence_vertices.erase(
            std::unique(parent.influence_vertices.begin(), parent.influence_vertices.end()),
            parent.influence_vertices.end());
        if (parent.influence_vertices.size() == num_graph_nodes) {
            parent.influence_is_full_graph = true;
            parent.influence_vertices.clear();
        }
    }
}

Mwpm& cluster_mwpm(ClusterSubgraph& cluster_subgraph, bool edge_correlations) {
    return cluster_subgraph_mwpm(cluster_subgraph, edge_correlations);
}

std::vector<uint8_t> build_child_watch_mask_in_parent(
    const ClusterSubgraph& parent_subgraph,
    const ClusterSubgraph& child_subgraph) {
    std::vector<uint8_t> watched(parent_subgraph.local_to_global_node_ids.size(), 0);
    for (size_t child_local_node_id = 0; child_local_node_id < child_subgraph.local_to_global_node_ids.size(); child_local_node_id++) {
        if (cluster_subgraph_local_node_is_boundary(child_subgraph, child_local_node_id)) {
            continue;
        }
        size_t global_node_id = child_subgraph.local_to_global_node_ids[child_local_node_id];
        if (global_node_id >= parent_subgraph.global_to_local_node_ids.size()) {
            continue;
        }
        size_t parent_local_node_id = parent_subgraph.global_to_local_node_ids[global_node_id];
        if (parent_local_node_id == SIZE_MAX || cluster_subgraph_local_node_is_boundary(parent_subgraph, parent_local_node_id)) {
            continue;
        }
        watched[parent_local_node_id] = 1;
    }
    return watched;
}

bool cluster_has_scheduled_work(ClusterSubgraph& cluster_subgraph, bool edge_correlations) {
    return cluster_mwpm(cluster_subgraph, edge_correlations).flooder.has_valid_tentative_events();
}

std::optional<cumulative_time_int> cluster_next_scheduled_time(
    ClusterSubgraph& cluster_subgraph, bool edge_correlations) {
    auto& flooder = cluster_mwpm(cluster_subgraph, edge_correlations).flooder;
    if (!flooder.has_valid_tentative_events()) {
        return std::nullopt;
    }
    try {
        return flooder.peek_next_valid_tentative_event_time();
    } catch (const std::invalid_argument&) {
        return std::nullopt;
    }
}

bool has_running_cluster_with_scheduled_work_at_time(
    LockstepDecodeResult& result, const LockstepSchedulerConfig& config, cumulative_time_int time) {
    for (size_t cluster_id = 0; cluster_id < result.cluster_subgraphs.size(); cluster_id++) {
        if (result.cluster_execution_states[cluster_id] != LockstepClusterExecutionState::RUNNING) {
            continue;
        }
        auto next_time = cluster_next_scheduled_time(result.cluster_subgraphs[cluster_id], config.edge_correlations);
        if (!next_time.has_value()) {
            continue;
        }
        if (*next_time == time) {
            return true;
        }
    }
    return false;
}

size_t checked_num_workers(const LockstepSchedulerConfig& config) {
    if (config.num_workers == 0) {
        throw std::invalid_argument("LockstepSchedulerConfig.num_workers must be at least 1.");
    }
    return config.num_workers;
}

bool has_running_cluster_with_scheduled_work(
    LockstepDecodeResult& result, const LockstepSchedulerConfig& config) {
    for (size_t cluster_id = 0; cluster_id < result.cluster_subgraphs.size(); cluster_id++) {
        if (result.cluster_execution_states[cluster_id] != LockstepClusterExecutionState::RUNNING) {
            continue;
        }
        if (cluster_has_scheduled_work(result.cluster_subgraphs[cluster_id], config.edge_correlations)) {
            return true;
        }
    }
    return false;
}

struct ReadyClusterScan {
    bool found = false;
    cumulative_time_int time = 0;
    std::vector<size_t> cluster_ids;
};

bool cluster_has_unimported_direct_child(
    const LockstepDecodeResult& result,
    size_t parent_cluster_id);

ReadyClusterScan scan_clusters_at_next_global_time(
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config) {
    ReadyClusterScan scan;
    for (size_t cluster_id = 0; cluster_id < result.cluster_subgraphs.size(); cluster_id++) {
        if (result.cluster_execution_states[cluster_id] != LockstepClusterExecutionState::RUNNING) {
            continue;
        }
        auto& cluster_subgraph = result.cluster_subgraphs[cluster_id];
        auto next_time = cluster_next_scheduled_time(cluster_subgraph, config.edge_correlations);
        if (!next_time.has_value()) {
            continue;
        }
        if (!scan.found || *next_time < scan.time) {
            scan.found = true;
            scan.time = *next_time;
            scan.cluster_ids.clear();
            scan.cluster_ids.push_back(cluster_id);
            continue;
        }
        if (*next_time == scan.time) {
            scan.cluster_ids.push_back(cluster_id);
        }
    }
    return scan;
}

bool maybe_mark_cluster_stopped(
    LockstepDecodeResult& result, size_t cluster_id, const LockstepSchedulerConfig& config) {
    auto previous_state = result.cluster_execution_states[cluster_id];
    bool previous_has_state = result.cluster_has_stopped_state[cluster_id];
    bool previous_provisional = result.cluster_result_is_provisional[cluster_id];
    bool previous_has_direct_result =
        (cluster_id < result.cluster_has_direct_result.size() && result.cluster_has_direct_result[cluster_id]) ||
        (cluster_id < result.cluster_has_direct_compact_result.size() &&
         result.cluster_has_direct_compact_result[cluster_id]);
    const bool is_root_cluster =
        result.hierarchy.clusters[cluster_id].parent_id == NO_PROCESSING_CLUSTER_PARENT;

    // Once a cluster's terminal export is already available, do not repeatedly
    // re-check and re-export it.  For roots this matters because the compact
    // matching extraction is destructive, and for all clusters it avoids charging
    // repeated stop polling to the critical-path lifecycle.  A root that is stopped
    // but still waiting for a direct child import intentionally falls through so it
    // can export after the last child arrives.
    if (previous_state == LockstepClusterExecutionState::QUIESCENT_STOPPED) {
        if ((!is_root_cluster && previous_has_state) || (is_root_cluster && previous_has_direct_result)) {
            return false;
        }
    }

    auto& subgraph = result.cluster_subgraphs[cluster_id];
    auto& mwpm = cluster_mwpm(subgraph, config.edge_correlations);
    if (mwpm.is_quiescent_without_shatter()) {
        result.cluster_execution_states[cluster_id] = LockstepClusterExecutionState::QUIESCENT_STOPPED;
        result.cluster_has_stopped_state[cluster_id] = true;
        if (is_root_cluster) {
            // Root clusters are never imported into a parent.  Do not build the
            // heavy MwpmStoppedState snapshot here; it saves detector ownership,
            // region ancestry, shell areas, etc. only for child-to-parent import.
            // For a root, the only required export is the final matching /
            // observable result.  Extract it once the root has no direct child
            // still waiting to be imported, because extraction shatters blossoms
            // and should not run before a later child import may resume the root.
            result.cluster_result_is_provisional[cluster_id] = false;
            const bool already_has_direct_result =
                (cluster_id < result.cluster_has_direct_result.size() && result.cluster_has_direct_result[cluster_id]) ||
                (cluster_id < result.cluster_has_direct_compact_result.size() &&
                 result.cluster_has_direct_compact_result[cluster_id]);
            if (!already_has_direct_result && !cluster_has_unimported_direct_child(result, cluster_id)) {
                if (mwpm.flooder.graph.num_observables <= sizeof(pm::obs_int) * 8) {
                    if (cluster_id >= result.cluster_direct_compact_results.size()) {
                        result.cluster_direct_compact_results.resize(cluster_id + 1);
                        result.cluster_has_direct_compact_result.resize(cluster_id + 1, false);
                    }
                    auto root_extract_start = steady_clock::now();
                    result.cluster_direct_compact_results[cluster_id] = extract_compact_result_from_current_mwpm_state(
                        mwpm, result.cluster_extraction_detection_events[cluster_id]);
                    auto root_extract_elapsed = elapsed_ns(root_extract_start);
                    result.profiling_stats.root_extraction_wall_ns += root_extract_elapsed;
                    if (cluster_id < result.cluster_mark_stop_wall_ns.size()) {
                        result.cluster_mark_stop_wall_ns[cluster_id] += root_extract_elapsed;
                    }
                    result.cluster_has_direct_compact_result[cluster_id] = true;
                    if (subgraph.use_shared_full_graph_mwpm) {
                        mwpm.reset();
                    }
                } else {
                    if (cluster_id >= result.cluster_direct_results.size()) {
                        result.cluster_direct_results.resize(cluster_id + 1);
                        result.cluster_has_direct_result.resize(cluster_id + 1, false);
                    }
                    auto root_extract_start = steady_clock::now();
                    result.cluster_direct_results[cluster_id] = extract_result_from_current_mwpm_state(
                        mwpm, result.cluster_extraction_detection_events[cluster_id]);
                    auto root_extract_elapsed = elapsed_ns(root_extract_start);
                    result.profiling_stats.root_extraction_wall_ns += root_extract_elapsed;
                    if (cluster_id < result.cluster_mark_stop_wall_ns.size()) {
                        result.cluster_mark_stop_wall_ns[cluster_id] += root_extract_elapsed;
                    }
                    result.cluster_has_direct_result[cluster_id] = true;
                    if (subgraph.use_shared_full_graph_mwpm) {
                        mwpm.reset();
                    }
                }
            }
        } else {
            auto stopped_export_start = steady_clock::now();
            if (config.enable_parentless_broadcast_import) {
                // Parentless/broadcast mode does not rely on a precomputed
                // parent influence subgraph to crop child state.  Instead, the
                // child exports only the detector ownership whose flooder source
                // is one of this cluster's active detectors.  Upper full-graph
                // workers can then import the filtered stopped state directly.
                result.cluster_stopped_states[cluster_id] = export_mwpm_stopped_state_for_active_detectors(
                    mwpm, subgraph.local_active_detectors);
            } else {
                result.cluster_stopped_states[cluster_id] = export_mwpm_stopped_state(mwpm);
            }
            if (cluster_id < result.cluster_mark_stop_wall_ns.size()) {
                result.cluster_mark_stop_wall_ns[cluster_id] += elapsed_ns(stopped_export_start);
            }
            result.cluster_result_is_provisional[cluster_id] = true;
        }
    } else if (!mwpm.flooder.has_valid_tentative_events()) {
        result.cluster_execution_states[cluster_id] = LockstepClusterExecutionState::DRAINED_WITH_LIVE_ALT_TREE;
        result.encountered_non_quiescent_cluster = true;
    }

    return previous_state != result.cluster_execution_states[cluster_id] ||
           previous_has_state != result.cluster_has_stopped_state[cluster_id] ||
           previous_provisional != result.cluster_result_is_provisional[cluster_id] ||
           previous_has_direct_result !=
               ((cluster_id < result.cluster_has_direct_result.size() && result.cluster_has_direct_result[cluster_id]) ||
                (cluster_id < result.cluster_has_direct_compact_result.size() &&
                 result.cluster_has_direct_compact_result[cluster_id]));
}

void ensure_extraction_seen_for_cluster(LockstepDecodeResult& result, size_t cluster_id) {
    if (cluster_id >= result.cluster_extraction_detection_event_seen.size()) {
        result.cluster_extraction_detection_event_seen.resize(cluster_id + 1);
    }
    auto& seen = result.cluster_extraction_detection_event_seen[cluster_id];
    const auto node_count = result.cluster_subgraphs[cluster_id].local_to_global_node_ids.size();
    if (seen.size() != node_count) {
        seen.assign(node_count, 0);
        for (auto node_id : result.cluster_extraction_detection_events[cluster_id]) {
            if (node_id < seen.size()) {
                seen[node_id] = 1;
            }
        }
    }
}

void rebuild_extraction_seen_for_cluster(LockstepDecodeResult& result, size_t cluster_id) {
    if (cluster_id >= result.cluster_extraction_detection_event_seen.size()) {
        result.cluster_extraction_detection_event_seen.resize(cluster_id + 1);
    }
    auto& seen = result.cluster_extraction_detection_event_seen[cluster_id];
    seen.assign(result.cluster_subgraphs[cluster_id].local_to_global_node_ids.size(), 0);
    for (auto node_id : result.cluster_extraction_detection_events[cluster_id]) {
        if (node_id < seen.size()) {
            seen[node_id] = 1;
        }
    }
}

void merge_child_extraction_detection_events_into_parent(
    LockstepDecodeResult& result,
    size_t parent_cluster_id,
    size_t child_cluster_id) {
    auto imported_sources = remap_source_detection_events_from_stopped_state(
        result.cluster_subgraphs[parent_cluster_id],
        result.cluster_subgraphs[child_cluster_id],
        result.cluster_stopped_states[child_cluster_id]);
    auto& parent_sources = result.cluster_extraction_detection_events[parent_cluster_id];
    ensure_extraction_seen_for_cluster(result, parent_cluster_id);
    auto& seen = result.cluster_extraction_detection_event_seen[parent_cluster_id];
    for (auto source : imported_sources) {
        if (source >= seen.size()) {
            throw std::invalid_argument("Imported extraction source is outside parent cluster subgraph.");
        }
        if (!seen[source]) {
            seen[source] = 1;
            parent_sources.push_back(source);
        }
    }
}

std::vector<size_t> direct_children_of_parent(
    const ProcessingClusterHierarchy& hierarchy,
    size_t parent_cluster_id) {
    std::vector<size_t> child_ids;
    for (size_t cluster_id = 0; cluster_id < hierarchy.clusters.size(); cluster_id++) {
        if (hierarchy.clusters[cluster_id].parent_id == parent_cluster_id) {
            child_ids.push_back(cluster_id);
        }
    }
    return child_ids;
}

bool cluster_has_unimported_direct_child(
    const LockstepDecodeResult& result,
    size_t parent_cluster_id) {
    for (auto child_cluster_id : result.direct_child_cluster_ids[parent_cluster_id]) {
        if (result.cluster_execution_states[child_cluster_id] != LockstepClusterExecutionState::IMPORTED_INTO_PARENT) {
            return true;
        }
    }
    return false;
}


bool cluster_is_empty_terminal_root(
    const LockstepDecodeResult& result,
    size_t cluster_id) {
    if (cluster_id >= result.hierarchy.clusters.size()) {
        return false;
    }
    const auto& cluster = result.hierarchy.clusters[cluster_id];
    return cluster.parent_id == NO_PROCESSING_CLUSTER_PARENT && cluster.active_detectors.empty();
}

MatchingResult extract_or_get_child_compact_result_for_terminal_root(
    LockstepDecodeResult& result,
    size_t child_cluster_id,
    const LockstepSchedulerConfig& config) {
    if (child_cluster_id < result.cluster_has_direct_compact_result.size() &&
        result.cluster_has_direct_compact_result[child_cluster_id]) {
        return result.cluster_direct_compact_results[child_cluster_id];
    }
    if (!result.cluster_has_stopped_state[child_cluster_id] ||
        result.cluster_execution_states[child_cluster_id] != LockstepClusterExecutionState::QUIESCENT_STOPPED) {
        throw std::invalid_argument("Terminal import requires a quiescent stopped child cluster.");
    }
    auto& child_subgraph = result.cluster_subgraphs[child_cluster_id];
    auto& child_mwpm = cluster_mwpm(child_subgraph, config.edge_correlations);
    auto child_extract_start = steady_clock::now();
    auto compact = extract_compact_result_from_current_mwpm_state(
        child_mwpm,
        result.cluster_extraction_detection_events[child_cluster_id]);
    auto elapsed = elapsed_ns(child_extract_start);
    result.profiling_stats.root_extraction_wall_ns += elapsed;
    if (child_cluster_id < result.cluster_mark_stop_wall_ns.size()) {
        result.cluster_mark_stop_wall_ns[child_cluster_id] += elapsed;
    }
    if (child_cluster_id >= result.cluster_direct_compact_results.size()) {
        result.cluster_direct_compact_results.resize(child_cluster_id + 1);
        result.cluster_has_direct_compact_result.resize(child_cluster_id + 1, false);
    }
    result.cluster_direct_compact_results[child_cluster_id] = compact;
    result.cluster_has_direct_compact_result[child_cluster_id] = true;
    if (child_subgraph.use_shared_full_graph_mwpm) {
        child_mwpm.reset();
    }
    return compact;
}

ExtendedMatchingResult extract_or_get_child_extended_result_for_terminal_root(
    LockstepDecodeResult& result,
    size_t child_cluster_id,
    const LockstepSchedulerConfig& config) {
    if (child_cluster_id < result.cluster_has_direct_result.size() &&
        result.cluster_has_direct_result[child_cluster_id]) {
        return result.cluster_direct_results[child_cluster_id];
    }
    if (!result.cluster_has_stopped_state[child_cluster_id] ||
        result.cluster_execution_states[child_cluster_id] != LockstepClusterExecutionState::QUIESCENT_STOPPED) {
        throw std::invalid_argument("Terminal import requires a quiescent stopped child cluster.");
    }
    auto& child_subgraph = result.cluster_subgraphs[child_cluster_id];
    auto& child_mwpm = cluster_mwpm(child_subgraph, config.edge_correlations);
    auto child_extract_start = steady_clock::now();
    auto extended = extract_result_from_current_mwpm_state(
        child_mwpm,
        result.cluster_extraction_detection_events[child_cluster_id]);
    auto elapsed = elapsed_ns(child_extract_start);
    result.profiling_stats.root_extraction_wall_ns += elapsed;
    if (child_cluster_id < result.cluster_mark_stop_wall_ns.size()) {
        result.cluster_mark_stop_wall_ns[child_cluster_id] += elapsed;
    }
    if (child_cluster_id >= result.cluster_direct_results.size()) {
        result.cluster_direct_results.resize(child_cluster_id + 1);
        result.cluster_has_direct_result.resize(child_cluster_id + 1, false);
    }
    result.cluster_direct_results[child_cluster_id] = extended;
    result.cluster_has_direct_result[child_cluster_id] = true;
    if (child_subgraph.use_shared_full_graph_mwpm) {
        child_mwpm.reset();
    }
    return extended;
}

void terminal_import_children_into_empty_root(
    LockstepDecodeResult& result,
    size_t parent_cluster_id,
    const std::vector<size_t>& child_ids,
    const LockstepSchedulerConfig& config) {
    const bool can_use_compact =
        cluster_mwpm(result.cluster_subgraphs[parent_cluster_id], config.edge_correlations)
            .flooder.graph.num_observables <= sizeof(pm::obs_int) * 8;

    if (can_use_compact) {
        if (parent_cluster_id >= result.cluster_direct_compact_results.size()) {
            result.cluster_direct_compact_results.resize(parent_cluster_id + 1);
            result.cluster_has_direct_compact_result.resize(parent_cluster_id + 1, false);
        }
        if (!result.cluster_has_direct_compact_result[parent_cluster_id]) {
            result.cluster_direct_compact_results[parent_cluster_id] = MatchingResult();
            result.cluster_has_direct_compact_result[parent_cluster_id] = true;
        }
        for (auto child_cluster_id : child_ids) {
            result.cluster_direct_compact_results[parent_cluster_id] +=
                extract_or_get_child_compact_result_for_terminal_root(result, child_cluster_id, config);
            result.profiling_stats.child_import_count++;
            result.cluster_execution_states[child_cluster_id] = LockstepClusterExecutionState::IMPORTED_INTO_PARENT;
            result.cluster_result_is_provisional[child_cluster_id] = false;
        }
    } else {
        if (parent_cluster_id >= result.cluster_direct_results.size()) {
            result.cluster_direct_results.resize(parent_cluster_id + 1);
            result.cluster_has_direct_result.resize(parent_cluster_id + 1, false);
        }
        if (!result.cluster_has_direct_result[parent_cluster_id]) {
            const size_t num_obs = cluster_mwpm(result.cluster_subgraphs[parent_cluster_id], config.edge_correlations)
                                       .flooder.graph.num_observables;
            result.cluster_direct_results[parent_cluster_id] = ExtendedMatchingResult(num_obs);
            result.cluster_has_direct_result[parent_cluster_id] = true;
        }
        for (auto child_cluster_id : child_ids) {
            result.cluster_direct_results[parent_cluster_id] +=
                extract_or_get_child_extended_result_for_terminal_root(result, child_cluster_id, config);
            result.profiling_stats.child_import_count++;
            result.cluster_execution_states[child_cluster_id] = LockstepClusterExecutionState::IMPORTED_INTO_PARENT;
            result.cluster_result_is_provisional[child_cluster_id] = false;
        }
    }

    result.cluster_execution_states[parent_cluster_id] = LockstepClusterExecutionState::QUIESCENT_STOPPED;
    result.cluster_has_stopped_state[parent_cluster_id] = true;
    result.cluster_result_is_provisional[parent_cluster_id] = false;
}

void capture_overlap_history_snapshot_for_cluster(
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config,
    size_t cluster_id) {
    if (!config.enable_overlap_checkpoints) {
        return;
    }
    if (cluster_id >= result.direct_child_cluster_ids.size() ||
        result.direct_child_cluster_ids[cluster_id].empty() ||
        !cluster_has_unimported_direct_child(result, cluster_id)) {
        return;
    }
    auto start = steady_clock::now();
    auto& subgraph = result.cluster_subgraphs[cluster_id];
    auto& mwpm = cluster_mwpm(subgraph, config.edge_correlations);
    auto snapshot = export_mwpm_live_snapshot(mwpm);
    auto time = mwpm.flooder.queue.cur_time;
    auto& history = result.cluster_overlap_snapshot_histories[cluster_id];
    if (!history.empty() && history.back().time == time) {
        history.back() = LockstepOverlapSnapshotHistoryEntry{
            time,
            std::move(snapshot),
            result.cluster_extraction_detection_events[cluster_id],
        };
    } else {
        history.push_back(LockstepOverlapSnapshotHistoryEntry{
            time,
            std::move(snapshot),
            result.cluster_extraction_detection_events[cluster_id],
        });
    }
    result.profiling_stats.overlap_checkpoint_capture_count++;
    auto capture_elapsed = elapsed_ns(start);
    result.profiling_stats.overlap_checkpoint_capture_wall_ns += capture_elapsed;
    if (cluster_id < result.cluster_checkpoint_capture_wall_ns.size()) {
        result.cluster_checkpoint_capture_wall_ns[cluster_id] += capture_elapsed;
    }
}

void maybe_capture_overlap_checkpoint_for_parent(
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config,
    size_t parent_cluster_id) {
    capture_overlap_history_snapshot_for_cluster(result, config, parent_cluster_id);
}

const LockstepOverlapSnapshotHistoryEntry* latest_overlap_snapshot_at_or_before(
    const LockstepDecodeResult& result,
    size_t cluster_id,
    cumulative_time_int time) {
    if (cluster_id >= result.cluster_overlap_snapshot_histories.size()) {
        return nullptr;
    }
    const auto& history = result.cluster_overlap_snapshot_histories[cluster_id];
    const LockstepOverlapSnapshotHistoryEntry* best = nullptr;
    for (const auto& entry : history) {
        if (entry.time <= time) {
            if (best == nullptr || entry.time >= best->time) {
                best = &entry;
            }
        }
    }
    return best;
}

void restore_parent_from_child_overlap_checkpoint(
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config,
    size_t parent_cluster_id,
    size_t child_cluster_id) {
    auto start = steady_clock::now();
    auto& parent_subgraph = result.cluster_subgraphs[parent_cluster_id];
    auto& parent_mwpm = cluster_mwpm(parent_subgraph, config.edge_correlations);
    const auto child_stop_time = result.cluster_stopped_states[child_cluster_id].algorithmic_time;
    const auto* checkpoint = latest_overlap_snapshot_at_or_before(result, parent_cluster_id, child_stop_time);
    if (checkpoint == nullptr) {
        // Backward-compatible fallback for old per-child checkpoints.  New runs
        // should normally use cluster_overlap_snapshot_histories instead.
        if (!result.cluster_has_overlap_checkpoint[child_cluster_id]) {
            throw std::invalid_argument(
                "Stopped-state import would overwrite detector ownership and no overlap checkpoint was captured.");
        }
        result.cluster_overlap_snapshot_histories[parent_cluster_id].push_back(LockstepOverlapSnapshotHistoryEntry{
            result.cluster_overlap_checkpoint_times[child_cluster_id],
            result.cluster_overlap_checkpoints[child_cluster_id],
            result.cluster_overlap_checkpoint_extraction_detection_events[child_cluster_id],
        });
        checkpoint = latest_overlap_snapshot_at_or_before(result, parent_cluster_id, child_stop_time);
        if (checkpoint == nullptr) {
            throw std::invalid_argument(
                "Stopped-state import would overwrite detector ownership and no compatible overlap checkpoint was captured.");
        }
    }
    const auto checkpoint_time = checkpoint->time;

    restore_mwpm_from_live_snapshot(parent_mwpm, checkpoint->snapshot);
    result.cluster_extraction_detection_events[parent_cluster_id] = checkpoint->extraction_detection_events;
    rebuild_extraction_seen_for_cluster(result, parent_cluster_id);

    std::vector<size_t> siblings_to_replay;
    for (auto sibling_cluster_id : result.direct_child_cluster_ids[parent_cluster_id]) {
        if (sibling_cluster_id == child_cluster_id) {
            continue;
        }
        if (result.cluster_execution_states[sibling_cluster_id] != LockstepClusterExecutionState::IMPORTED_INTO_PARENT ||
            !result.cluster_has_stopped_state[sibling_cluster_id]) {
            continue;
        }
        const auto sibling_stop_time = result.cluster_stopped_states[sibling_cluster_id].algorithmic_time;
        if (sibling_stop_time <= checkpoint_time || sibling_stop_time > child_stop_time) {
            continue;
        }
        siblings_to_replay.push_back(sibling_cluster_id);
    }
    std::sort(siblings_to_replay.begin(), siblings_to_replay.end(), [&](size_t a, size_t b) {
        auto ta = result.cluster_stopped_states[a].algorithmic_time;
        auto tb = result.cluster_stopped_states[b].algorithmic_time;
        if (ta != tb) return ta < tb;
        return a < b;
    });
    for (auto sibling_cluster_id : siblings_to_replay) {
        import_stopped_state_into_cluster_subgraph(
            parent_subgraph,
            result.cluster_subgraphs[sibling_cluster_id],
            result.cluster_stopped_states[sibling_cluster_id],
            config.edge_correlations);
        merge_child_extraction_detection_events_into_parent(result, parent_cluster_id, sibling_cluster_id);
    }
    result.cluster_used_overlap_checkpoint_recovery[child_cluster_id] = true;
    result.profiling_stats.overlap_checkpoint_restore_count++;
    auto restore_elapsed = elapsed_ns(start);
    result.profiling_stats.overlap_checkpoint_restore_wall_ns += restore_elapsed;
    if (parent_cluster_id < result.cluster_checkpoint_restore_wall_ns.size()) {
        result.cluster_checkpoint_restore_wall_ns[parent_cluster_id] += restore_elapsed;
    }
}

void restart_parent_from_time_zero_and_replay_imported_siblings(
    UserGraph& graph,
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config,
    size_t parent_cluster_id,
    size_t child_cluster_id) {
    auto start = steady_clock::now();
    const auto child_stop_time = result.cluster_stopped_states[child_cluster_id].algorithmic_time;

    // Fast-path recovery when checkpoint capture is disabled:
    // rebuild the parent cluster MWPM from algorithmic time 0 and replay only
    // children that had already been imported before the new child.  This keeps
    // the heavyweight checkpoint implementation available for experiments while
    // making the default runtime path pay zero checkpoint-capture cost.
    result.cluster_subgraphs[parent_cluster_id] =
        build_worker_cluster_state(graph, result.hierarchy.clusters[parent_cluster_id], config);
    auto& parent_subgraph = result.cluster_subgraphs[parent_cluster_id];
    auto& parent_mwpm = cluster_mwpm(parent_subgraph, config.edge_correlations);
    result.cluster_extraction_detection_events[parent_cluster_id] = parent_subgraph.local_active_detectors;
    rebuild_extraction_seen_for_cluster(result, parent_cluster_id);
    initialize_mwpm_for_detection_events(parent_mwpm, parent_subgraph.local_active_detectors);
    result.cluster_execution_states[parent_cluster_id] = LockstepClusterExecutionState::RUNNING;
    result.cluster_has_stopped_state[parent_cluster_id] = false;
    result.cluster_result_is_provisional[parent_cluster_id] = false;
    if (parent_cluster_id < result.cluster_has_direct_result.size()) {
        result.cluster_has_direct_result[parent_cluster_id] = false;
    }
    if (parent_cluster_id < result.cluster_overlap_snapshot_histories.size()) {
        result.cluster_overlap_snapshot_histories[parent_cluster_id].clear();
    }

    std::vector<size_t> siblings_to_replay;
    for (auto sibling_cluster_id : result.direct_child_cluster_ids[parent_cluster_id]) {
        if (sibling_cluster_id == child_cluster_id) {
            continue;
        }
        if (result.cluster_execution_states[sibling_cluster_id] != LockstepClusterExecutionState::IMPORTED_INTO_PARENT ||
            !result.cluster_has_stopped_state[sibling_cluster_id]) {
            continue;
        }
        const auto sibling_stop_time = result.cluster_stopped_states[sibling_cluster_id].algorithmic_time;
        if (sibling_stop_time > child_stop_time) {
            continue;
        }
        siblings_to_replay.push_back(sibling_cluster_id);
    }
    std::sort(siblings_to_replay.begin(), siblings_to_replay.end(), [&](size_t a, size_t b) {
        auto ta = result.cluster_stopped_states[a].algorithmic_time;
        auto tb = result.cluster_stopped_states[b].algorithmic_time;
        if (ta != tb) return ta < tb;
        return a < b;
    });
    for (auto sibling_cluster_id : siblings_to_replay) {
        import_stopped_state_into_cluster_subgraph(
            parent_subgraph,
            result.cluster_subgraphs[sibling_cluster_id],
            result.cluster_stopped_states[sibling_cluster_id],
            config.edge_correlations);
        merge_child_extraction_detection_events_into_parent(result, parent_cluster_id, sibling_cluster_id);
    }

    result.cluster_used_overlap_checkpoint_recovery[child_cluster_id] = true;
    result.profiling_stats.overlap_checkpoint_restore_count++;
    auto restart_elapsed = elapsed_ns(start);
    result.profiling_stats.overlap_checkpoint_restore_wall_ns += restart_elapsed;
    if (parent_cluster_id < result.cluster_checkpoint_restore_wall_ns.size()) {
        result.cluster_checkpoint_restore_wall_ns[parent_cluster_id] += restart_elapsed;
    }
}


void debug_print_import_overwrite_details(
    const LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config,
    size_t parent_cluster_id,
    size_t child_cluster_id,
    const char* phase) {
    const auto& target_subgraph = result.cluster_subgraphs[parent_cluster_id];
    const auto& source_subgraph = result.cluster_subgraphs[child_cluster_id];
    const auto& source_state = result.cluster_stopped_states[child_cluster_id];
    auto& target_mwpm = cluster_mwpm(const_cast<ClusterSubgraph&>(target_subgraph), config.edge_correlations);
    size_t printed = 0;
    size_t conflicts = 0;
    for (const auto& detector_state : source_state.detector_states) {
        if (detector_state.node_id == NO_STOPPED_NODE || detector_state.region_id == NO_STOPPED_REGION) {
            continue;
        }
        if (detector_state.node_id >= source_subgraph.local_to_global_node_ids.size()) {
            continue;
        }
        size_t global_node_id = source_subgraph.local_to_global_node_ids[detector_state.node_id];
        if (global_node_id >= target_subgraph.global_to_local_node_ids.size()) {
            continue;
        }
        size_t target_local_node_id = target_subgraph.global_to_local_node_ids[global_node_id];
        if (target_local_node_id == SIZE_MAX || target_local_node_id >= target_mwpm.flooder.graph.nodes.size()) {
            continue;
        }
        const auto& target_node = target_mwpm.flooder.graph.nodes[target_local_node_id];
        if (target_node.region_that_arrived == nullptr) {
            continue;
        }
        conflicts++;
        if (printed < 12) {
            bool target_is_parent_source = std::find(
                target_subgraph.local_active_detectors.begin(),
                target_subgraph.local_active_detectors.end(),
                (uint64_t)target_local_node_id) != target_subgraph.local_active_detectors.end();
            bool source_is_child_source = detector_state.reached_from_source_id == detector_state.node_id;
            const auto* chosen_checkpoint = latest_overlap_snapshot_at_or_before(
                result, parent_cluster_id, source_state.algorithmic_time);
            auto chosen_checkpoint_time = chosen_checkpoint == nullptr ? (cumulative_time_int)0 : chosen_checkpoint->time;
            std::cerr << "IMPORT_CONFLICT phase=" << phase
                      << " parent=" << parent_cluster_id
                      << " child=" << child_cluster_id
                      << " parent_level=" << result.hierarchy.clusters[parent_cluster_id].level
                      << " child_level=" << result.hierarchy.clusters[child_cluster_id].level
                      << " parent_state=" << (int)result.cluster_execution_states[parent_cluster_id]
                      << " child_state=" << (int)result.cluster_execution_states[child_cluster_id]
                      << " child_stop_time=" << source_state.algorithmic_time
                      << " chosen_checkpoint_time=" << chosen_checkpoint_time
                      << " parent_cur_time=" << target_mwpm.flooder.queue.cur_time
                      << " global_node=" << global_node_id
                      << " target_local=" << target_local_node_id
                      << " source_local=" << detector_state.node_id
                      << " target_is_parent_source=" << target_is_parent_source
                      << " source_is_child_source=" << source_is_child_source
                      << " target_region=" << (void*)target_node.region_that_arrived
                      << " target_reached_from=" << (target_node.reached_from_source == nullptr ? -1 : (long long)(target_node.reached_from_source - target_mwpm.flooder.graph.nodes.data()))
                      << "\n";
            printed++;
        }
    }
    if (conflicts > 0) {
        std::cerr << "IMPORT_CONFLICT_SUMMARY phase=" << phase
                  << " parent=" << parent_cluster_id
                  << " child=" << child_cluster_id
                  << " conflicts=" << conflicts << "\n";
    }
}

bool import_newly_stopped_children_into_parents_at_time(
    UserGraph& graph,
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config) {
    bool imported_any_ever = false;
    bool imported_any = true;
    while (imported_any) {
        imported_any = false;
        for (size_t child_cluster_id = 0; child_cluster_id < result.cluster_subgraphs.size(); child_cluster_id++) {
            if (result.cluster_execution_states[child_cluster_id] != LockstepClusterExecutionState::QUIESCENT_STOPPED) {
                continue;
            }
            auto parent_cluster_id = result.hierarchy.clusters[child_cluster_id].parent_id;
            if (parent_cluster_id == NO_PROCESSING_CLUSTER_PARENT) {
                continue;
            }
            auto& parent_subgraph = result.cluster_subgraphs[parent_cluster_id];
            auto& parent_mwpm = cluster_mwpm(parent_subgraph, config.edge_correlations);
            if (stopped_state_import_would_overwrite_detector_ownership(
                    parent_subgraph,
                    result.cluster_subgraphs[child_cluster_id],
                    result.cluster_stopped_states[child_cluster_id])) {
                debug_print_import_overwrite_details(result, config, parent_cluster_id, child_cluster_id, "before_recovery");
                if (config.enable_overlap_checkpoints) {
                    restore_parent_from_child_overlap_checkpoint(result, config, parent_cluster_id, child_cluster_id);
                } else {
                    restart_parent_from_time_zero_and_replay_imported_siblings(
                        graph, result, config, parent_cluster_id, child_cluster_id);
                }
                if (stopped_state_import_would_overwrite_detector_ownership(
                        parent_subgraph,
                        result.cluster_subgraphs[child_cluster_id],
                        result.cluster_stopped_states[child_cluster_id])) {
                    debug_print_import_overwrite_details(result, config, parent_cluster_id, child_cluster_id, "after_recovery");
                }
            } else {
                parent_mwpm.flooder.queue.cur_time = std::max(
                    parent_mwpm.flooder.queue.cur_time,
                    result.cluster_stopped_states[child_cluster_id].algorithmic_time);
            }

            auto import_start = steady_clock::now();
            import_stopped_state_into_cluster_subgraph(
                parent_subgraph,
                result.cluster_subgraphs[child_cluster_id],
                result.cluster_stopped_states[child_cluster_id],
                config.edge_correlations);
            merge_child_extraction_detection_events_into_parent(result, parent_cluster_id, child_cluster_id);
            result.profiling_stats.child_import_count++;
            auto import_elapsed = elapsed_ns(import_start);
            result.profiling_stats.child_import_wall_ns += import_elapsed;
            if (parent_cluster_id < result.cluster_child_import_wall_ns.size()) {
                result.cluster_child_import_wall_ns[parent_cluster_id] += import_elapsed;
            }
            result.cluster_execution_states[child_cluster_id] = LockstepClusterExecutionState::IMPORTED_INTO_PARENT;
            result.cluster_result_is_provisional[child_cluster_id] = false;

            if (result.cluster_execution_states[parent_cluster_id] == LockstepClusterExecutionState::QUIESCENT_STOPPED) {
                result.cluster_execution_states[parent_cluster_id] = LockstepClusterExecutionState::RUNNING;
                result.cluster_has_stopped_state[parent_cluster_id] = false;
            } else if (result.cluster_execution_states[parent_cluster_id] == LockstepClusterExecutionState::DRAINED_WITH_LIVE_ALT_TREE) {
                result.cluster_execution_states[parent_cluster_id] = LockstepClusterExecutionState::RUNNING;
            }
            maybe_mark_cluster_stopped(result, parent_cluster_id, config);
            capture_overlap_history_snapshot_for_cluster(result, config, parent_cluster_id);
            imported_any = true;
            imported_any_ever = true;
        }
    }
    return imported_any_ever;
}


bool level_is_ready_for_batched_child_import(
    const LockstepDecodeResult& result,
    size_t child_level) {
    bool any_waiting = false;
    for (size_t child_cluster_id = 0; child_cluster_id < result.hierarchy.clusters.size(); child_cluster_id++) {
        const auto& child = result.hierarchy.clusters[child_cluster_id];
        if (child.level != child_level || child.parent_id == NO_PROCESSING_CLUSTER_PARENT) {
            continue;
        }
        auto state = result.cluster_execution_states[child_cluster_id];
        if (state == LockstepClusterExecutionState::QUIESCENT_STOPPED) {
            any_waiting = true;
            continue;
        }
        if (state == LockstepClusterExecutionState::IMPORTED_INTO_PARENT) {
            continue;
        }
        return false;
    }
    return any_waiting;
}

bool import_ready_level_batched_children_into_parents(
    UserGraph& graph,
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config) {
    (void)graph;
    bool imported_any_ever = false;
    const size_t max_level = find_maximum_level(result.hierarchy);
    for (size_t child_level = 0; child_level < max_level; child_level++) {
        if (!level_is_ready_for_batched_child_import(result, child_level)) {
            continue;
        }

        std::vector<std::vector<size_t>> children_by_parent(result.hierarchy.clusters.size());
        for (size_t child_cluster_id = 0; child_cluster_id < result.hierarchy.clusters.size(); child_cluster_id++) {
            const auto& child = result.hierarchy.clusters[child_cluster_id];
            if (child.level != child_level || child.parent_id == NO_PROCESSING_CLUSTER_PARENT) {
                continue;
            }
            if (result.cluster_execution_states[child_cluster_id] == LockstepClusterExecutionState::QUIESCENT_STOPPED) {
                children_by_parent[child.parent_id].push_back(child_cluster_id);
            }
        }

        for (size_t parent_cluster_id = 0; parent_cluster_id < children_by_parent.size(); parent_cluster_id++) {
            auto& child_ids = children_by_parent[parent_cluster_id];
            if (child_ids.empty()) {
                continue;
            }
            std::sort(child_ids.begin(), child_ids.end(), [&](size_t a, size_t b) {
                auto ta = result.cluster_stopped_states[a].algorithmic_time;
                auto tb = result.cluster_stopped_states[b].algorithmic_time;
                if (ta != tb) return ta < tb;
                return a < b;
            });
            auto& parent_subgraph = result.cluster_subgraphs[parent_cluster_id];
            auto import_start = steady_clock::now();
            if (cluster_is_empty_terminal_root(result, parent_cluster_id)) {
                // Terminal aggregation fast path.  A forced single final root with
                // no active detectors does not need the child flooder states: no
                // later parent growth can change the already-stopped child
                // matchings.  Aggregate the child matching/observable results and
                // skip stopped-state import, detector ownership replay, and
                // rescheduling entirely.
                terminal_import_children_into_empty_root(result, parent_cluster_id, child_ids, config);
            } else {
                std::vector<size_t> imported_target_local_node_ids;
                for (auto child_cluster_id : child_ids) {
                    StoppedStateImportOptions options;
                    options.reschedule_after_import = false;
                    options.imported_target_local_node_ids = &imported_target_local_node_ids;
                    import_stopped_state_into_cluster_subgraph(
                        parent_subgraph,
                        result.cluster_subgraphs[child_cluster_id],
                        result.cluster_stopped_states[child_cluster_id],
                        options,
                        config.edge_correlations);
                    merge_child_extraction_detection_events_into_parent(result, parent_cluster_id, child_cluster_id);
                    result.profiling_stats.child_import_count++;
                    result.cluster_execution_states[child_cluster_id] = LockstepClusterExecutionState::IMPORTED_INTO_PARENT;
                    result.cluster_result_is_provisional[child_cluster_id] = false;
                }
                // Reschedule once for this parent after all level-k children that target it have been imported.
                // In particular, a root receives at most one reschedule batch per lower level, hence at most k_max-1.
                reschedule_imported_detector_nodes(parent_subgraph, imported_target_local_node_ids, config.edge_correlations);
                if (result.cluster_execution_states[parent_cluster_id] == LockstepClusterExecutionState::QUIESCENT_STOPPED) {
                    result.cluster_execution_states[parent_cluster_id] = LockstepClusterExecutionState::RUNNING;
                    result.cluster_has_stopped_state[parent_cluster_id] = false;
                } else if (result.cluster_execution_states[parent_cluster_id] == LockstepClusterExecutionState::DRAINED_WITH_LIVE_ALT_TREE) {
                    result.cluster_execution_states[parent_cluster_id] = LockstepClusterExecutionState::RUNNING;
                }
                maybe_mark_cluster_stopped(result, parent_cluster_id, config);
                capture_overlap_history_snapshot_for_cluster(result, config, parent_cluster_id);
            }
            auto import_elapsed = elapsed_ns(import_start);
            result.profiling_stats.child_import_wall_ns += import_elapsed;
            if (parent_cluster_id < result.cluster_child_import_wall_ns.size()) {
                result.cluster_child_import_wall_ns[parent_cluster_id] += import_elapsed;
            }
            imported_any_ever = true;
        }
    }
    return imported_any_ever;
}


bool level_is_ready_for_parentless_broadcast_import(
    const LockstepDecodeResult& result,
    size_t child_level) {
    bool any_waiting = false;
    bool any_upper = false;
    for (const auto& cluster : result.hierarchy.clusters) {
        if (cluster.level > child_level) {
            any_upper = true;
            break;
        }
    }
    if (!any_upper) {
        return false;
    }
    for (size_t child_cluster_id = 0; child_cluster_id < result.hierarchy.clusters.size(); child_cluster_id++) {
        const auto& child = result.hierarchy.clusters[child_cluster_id];
        if (child.level != child_level) {
            continue;
        }
        auto state = result.cluster_execution_states[child_cluster_id];
        if (state == LockstepClusterExecutionState::QUIESCENT_STOPPED) {
            any_waiting = true;
            continue;
        }
        if (state == LockstepClusterExecutionState::IMPORTED_INTO_PARENT) {
            continue;
        }
        return false;
    }
    return any_waiting;
}

std::vector<size_t> parentless_broadcast_targets_for_level(
    const LockstepDecodeResult& result,
    size_t child_level) {
    std::vector<size_t> targets;
    for (size_t cluster_id = 0; cluster_id < result.hierarchy.clusters.size(); cluster_id++) {
        const auto& cluster = result.hierarchy.clusters[cluster_id];
        if (cluster.level <= child_level) {
            continue;
        }
        if (!cluster.active_detectors.empty()) {
            targets.push_back(cluster_id);
        }
    }
    if (!targets.empty()) {
        return targets;
    }
    // If there is no non-empty upper cluster, deliver to terminal empty roots.
    // This preserves the existing final aggregation step without sending every
    // lower level directly to the terminal root, which would double count states
    // already inherited through an intermediate upper cluster.
    for (size_t cluster_id = 0; cluster_id < result.hierarchy.clusters.size(); cluster_id++) {
        const auto& cluster = result.hierarchy.clusters[cluster_id];
        if (cluster.level > child_level && cluster.active_detectors.empty()) {
            targets.push_back(cluster_id);
        }
    }
    return targets;
}

bool import_ready_level_parentless_broadcast(
    UserGraph& graph,
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config) {
    (void)graph;
    bool imported_any_ever = false;
    const size_t max_level = find_maximum_level(result.hierarchy);
    for (size_t child_level = 0; child_level < max_level; child_level++) {
        if (!level_is_ready_for_parentless_broadcast_import(result, child_level)) {
            continue;
        }
        std::vector<size_t> child_ids;
        for (size_t child_cluster_id = 0; child_cluster_id < result.hierarchy.clusters.size(); child_cluster_id++) {
            const auto& child = result.hierarchy.clusters[child_cluster_id];
            if (child.level != child_level) {
                continue;
            }
            if (result.cluster_execution_states[child_cluster_id] == LockstepClusterExecutionState::QUIESCENT_STOPPED) {
                child_ids.push_back(child_cluster_id);
            }
        }
        if (child_ids.empty()) {
            continue;
        }
        std::sort(child_ids.begin(), child_ids.end(), [&](size_t a, size_t b) {
            auto ta = result.cluster_stopped_states[a].algorithmic_time;
            auto tb = result.cluster_stopped_states[b].algorithmic_time;
            if (ta != tb) return ta < tb;
            return a < b;
        });

        auto targets = parentless_broadcast_targets_for_level(result, child_level);
        for (auto target_cluster_id : targets) {
            auto& target_subgraph = result.cluster_subgraphs[target_cluster_id];
            auto import_start = steady_clock::now();
            if (cluster_is_empty_terminal_root(result, target_cluster_id)) {
                terminal_import_children_into_empty_root(result, target_cluster_id, child_ids, config);
            } else {
                std::vector<size_t> imported_target_local_node_ids;
                for (auto child_cluster_id : child_ids) {
                    StoppedStateImportOptions options;
                    options.reschedule_after_import = false;
                    options.imported_target_local_node_ids = &imported_target_local_node_ids;
                    // Parentless full-graph workers rely on export-side active-detector
                    // filtering.  Import-side filtering would require a target influence
                    // subgraph, which the new clustering path intentionally does not
                    // build.  Duplicate arrivals from multiple upper-level routes are
                    // allowed to overwrite; theoretically harmful overwrites should not
                    // occur for the active-detector-limited export.
                    options.skip_unmappable_or_conflicting_state = false;
                    options.allow_detector_ownership_overwrite = true;
                    import_stopped_state_into_cluster_subgraph(
                        target_subgraph,
                        result.cluster_subgraphs[child_cluster_id],
                        result.cluster_stopped_states[child_cluster_id],
                        options,
                        config.edge_correlations);
                    if (!imported_target_local_node_ids.empty()) {
                        merge_child_extraction_detection_events_into_parent(result, target_cluster_id, child_cluster_id);
                    }
                    result.profiling_stats.child_import_count++;
                }
                reschedule_imported_detector_nodes(target_subgraph, imported_target_local_node_ids, config.edge_correlations);
                if (!imported_target_local_node_ids.empty()) {
                    if (result.cluster_execution_states[target_cluster_id] == LockstepClusterExecutionState::QUIESCENT_STOPPED) {
                        result.cluster_execution_states[target_cluster_id] = LockstepClusterExecutionState::RUNNING;
                        result.cluster_has_stopped_state[target_cluster_id] = false;
                    } else if (result.cluster_execution_states[target_cluster_id] == LockstepClusterExecutionState::DRAINED_WITH_LIVE_ALT_TREE) {
                        result.cluster_execution_states[target_cluster_id] = LockstepClusterExecutionState::RUNNING;
                    }
                    maybe_mark_cluster_stopped(result, target_cluster_id, config);
                    capture_overlap_history_snapshot_for_cluster(result, config, target_cluster_id);
                }
            }
            auto import_elapsed = elapsed_ns(import_start);
            result.profiling_stats.child_import_wall_ns += import_elapsed;
            if (target_cluster_id < result.cluster_child_import_wall_ns.size()) {
                result.cluster_child_import_wall_ns[target_cluster_id] += import_elapsed;
            }
            imported_any_ever = true;
        }
        for (auto child_cluster_id : child_ids) {
            result.cluster_execution_states[child_cluster_id] = LockstepClusterExecutionState::IMPORTED_INTO_PARENT;
            result.cluster_result_is_provisional[child_cluster_id] = false;
        }
    }
    return imported_any_ever;
}

bool mark_stopped_running_clusters_and_import_until_stable(
    UserGraph& graph,
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config) {
    bool progressed_ever = false;
    while (true) {
        bool progressed = false;
        for (size_t cluster_id = 0; cluster_id < result.cluster_subgraphs.size(); cluster_id++) {
            if (result.cluster_execution_states[cluster_id] != LockstepClusterExecutionState::RUNNING) {
                continue;
            }
            auto& mwpm = cluster_mwpm(result.cluster_subgraphs[cluster_id], config.edge_correlations);
            if (!mwpm.flooder.has_valid_tentative_events()) {
                progressed |= maybe_mark_cluster_stopped(result, cluster_id, config);
            }
        }
        if (config.enable_parentless_broadcast_import) {
            progressed |= import_ready_level_parentless_broadcast(graph, result, config);
        } else if (config.enable_level_batched_scheduler) {
            progressed |= import_ready_level_batched_children_into_parents(graph, result, config);
        } else {
            progressed |= import_newly_stopped_children_into_parents_at_time(graph, result, config);
        }
        progressed_ever |= progressed;
        if (!progressed) {
            break;
        }
    }
    return progressed_ever;
}

struct ClusterStepOutcome {
    size_t cluster_id;
    cumulative_time_int progress_time;
};

std::vector<size_t> collect_clusters_to_step_at_time(
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config,
    cumulative_time_int next_global_time) {
    std::vector<size_t> cluster_ids;
    for (size_t cluster_id = 0; cluster_id < result.cluster_subgraphs.size(); cluster_id++) {
        if (result.cluster_execution_states[cluster_id] != LockstepClusterExecutionState::RUNNING) {
            continue;
        }
        auto& mwpm = cluster_mwpm(result.cluster_subgraphs[cluster_id], config.edge_correlations);
        auto next_time = cluster_next_scheduled_time(result.cluster_subgraphs[cluster_id], config.edge_correlations);
        if (!next_time.has_value()) {
            continue;
        }
        if (*next_time == next_global_time) {
            cluster_ids.push_back(cluster_id);
        }
    }
    return cluster_ids;
}

ClusterStepOutcome step_cluster_through_time_slice(
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config,
    size_t cluster_id,
    cumulative_time_int time_slice) {
    auto& mwpm = cluster_mwpm(result.cluster_subgraphs[cluster_id], config.edge_correlations);
    auto next_time = cluster_next_scheduled_time(result.cluster_subgraphs[cluster_id], config.edge_correlations);
    if (next_time.has_value() && *next_time == time_slice) {
        auto cluster_step_start = steady_clock::now();
        try {
            auto bucket_event_count = flooder_queue_bucket_event_count(mwpm);
            if (mwpm.flooder.queue.empty() || bucket_event_count != mwpm.flooder.queue.size()) {
                auto old_qsize = mwpm.flooder.queue.size();
                auto requeued = rebuild_flooder_queue_from_desired_trackers(mwpm);
                std::cerr << "QUEUE_REBUILD cluster=" << cluster_id
                          << " cur=" << mwpm.flooder.queue.cur_time
                          << " old_qsize=" << old_qsize
                          << " old_bucket_count=" << bucket_event_count
                          << " requeued=" << requeued
                          << " qsize=" << mwpm.flooder.queue.size() << "\n";
            }
            auto event = mwpm.flooder.step_one_tentative_event_returning_mwpm_event();
            if (event.event_type != NO_EVENT) {
                mwpm.process_event(event);
            }
            if (cluster_id < result.cluster_step_wall_ns.size()) {
                result.cluster_step_wall_ns[cluster_id] += elapsed_ns(cluster_step_start);
            }
        } catch (const std::exception& ex) {
            if (cluster_id < result.cluster_step_wall_ns.size()) {
                result.cluster_step_wall_ns[cluster_id] += elapsed_ns(cluster_step_start);
            }
            const auto& cluster = result.hierarchy.clusters[cluster_id];
            throw std::invalid_argument(
                std::string("Failed while stepping cluster=") + std::to_string(cluster_id) +
                " level=" + std::to_string(cluster.level) +
                " parent=" + std::to_string(cluster.parent_id) +
                " active=" + std::to_string(cluster.active_detectors.size()) +
                " time_slice=" + std::to_string(time_slice) +
                " queue_cur=" + std::to_string(mwpm.flooder.queue.cur_time) +
                " queue_size=" + std::to_string(mwpm.flooder.queue.size()) +
                " error=" + ex.what());
        }
    }
    return {cluster_id, mwpm.flooder.queue.cur_time};
}

std::vector<ClusterStepOutcome> step_clusters_at_time(
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config,
    const std::vector<size_t>& cluster_ids_to_step,
    cumulative_time_int time_slice) {
    std::vector<ClusterStepOutcome> outcomes;
    if (cluster_ids_to_step.empty()) {
        return outcomes;
    }

    size_t num_workers = std::min(checked_num_workers(config), cluster_ids_to_step.size());
    if (num_workers <= 1 || cluster_ids_to_step.size() <= 1) {
        outcomes.reserve(cluster_ids_to_step.size());
        for (auto cluster_id : cluster_ids_to_step) {
            outcomes.push_back(step_cluster_through_time_slice(result, config, cluster_id, time_slice));
        }
        return outcomes;
    }

    std::vector<std::vector<ClusterStepOutcome>> worker_outcomes(num_workers);
    std::vector<std::future<void>> futures;
    futures.reserve(num_workers);
    size_t chunk_size = (cluster_ids_to_step.size() + num_workers - 1) / num_workers;
    for (size_t worker = 0; worker < num_workers; worker++) {
        size_t start = worker * chunk_size;
        if (start >= cluster_ids_to_step.size()) {
            break;
        }
        size_t end = std::min(cluster_ids_to_step.size(), start + chunk_size);
        futures.push_back(std::async(std::launch::async, [&, worker, start, end]() {
            auto& local_outcomes = worker_outcomes[worker];
            local_outcomes.reserve(end - start);
            for (size_t k = start; k < end; k++) {
                auto cluster_id = cluster_ids_to_step[k];
                local_outcomes.push_back(step_cluster_through_time_slice(result, config, cluster_id, time_slice));
            }
        }));
    }
    for (auto& future : futures) {
        future.get();
    }
    for (auto& local_outcomes : worker_outcomes) {
        outcomes.insert(outcomes.end(), local_outcomes.begin(), local_outcomes.end());
    }
    return outcomes;
}

size_t preprocessing_worker_count(size_t item_count) {
    if (item_count <= 1) {
        return 1;
    }
    size_t workers = 1;
    if (const char* env = std::getenv("PYMATCHING_PARALLEL_PREPROCESS_WORKERS")) {
        workers = std::max<size_t>(1, static_cast<size_t>(std::stoull(env)));
    }
    return std::min(workers, item_count);
}


ClusterSubgraph build_worker_cluster_state(
    UserGraph& graph,
    const ProcessingCluster& cluster,
    const LockstepSchedulerConfig& config,
    bool allow_shared_global_mwpm) {
    if (config.enable_full_graph_worker_subgraphs) {
        const bool lazy_acquire_persistent_mwpm =
            config.enable_parentless_level_pipeline &&
            config.enable_parentless_broadcast_import &&
            config.enable_growing_only_processing_clusters &&
            !config.edge_correlations;
        return build_full_graph_worker_cluster_subgraph(
            graph, cluster, allow_shared_global_mwpm, lazy_acquire_persistent_mwpm);
    }
    return build_cluster_subgraph(graph, cluster);
}

std::vector<ClusterSubgraph> build_cluster_subgraphs_parallel(
    UserGraph& graph,
    const std::vector<ProcessingCluster>& clusters,
    const LockstepSchedulerConfig& config,
    std::vector<uint64_t>* cluster_build_wall_ns = nullptr) {
    if (cluster_build_wall_ns != nullptr) {
        cluster_build_wall_ns->assign(clusters.size(), 0);
    }
    size_t shareable_root_count = 0;
    for (const auto& cluster : clusters) {
        if (cluster.parent_id == NO_PROCESSING_CLUSTER_PARENT && !cluster.forced_by_max_level) {
            shareable_root_count++;
        }
    }
    const bool allow_shared_global_mwpm = shareable_root_count == 1;

    const size_t worker_count = preprocessing_worker_count(clusters.size());
    if (worker_count <= 1) {
        std::vector<ClusterSubgraph> subgraphs;
        subgraphs.reserve(clusters.size());
        for (size_t k = 0; k < clusters.size(); k++) {
            auto build_start = steady_clock::now();
            subgraphs.push_back(build_worker_cluster_state(graph, clusters[k], config, allow_shared_global_mwpm));
            if (cluster_build_wall_ns != nullptr) {
                (*cluster_build_wall_ns)[k] = elapsed_ns(build_start);
            }
        }
        return subgraphs;
    }

    std::vector<std::optional<ClusterSubgraph>> slots(clusters.size());
    std::vector<std::future<void>> futures;
    futures.reserve(worker_count);
    const size_t chunk_size = (clusters.size() + worker_count - 1) / worker_count;
    for (size_t worker = 0; worker < worker_count; worker++) {
        const size_t start = worker * chunk_size;
        if (start >= clusters.size()) {
            break;
        }
        const size_t end = std::min(clusters.size(), start + chunk_size);
        futures.push_back(std::async(std::launch::async, [&, start, end]() {
            for (size_t k = start; k < end; k++) {
                auto build_start = steady_clock::now();
                slots[k].emplace(build_worker_cluster_state(graph, clusters[k], config, allow_shared_global_mwpm));
                if (cluster_build_wall_ns != nullptr) {
                    (*cluster_build_wall_ns)[k] = elapsed_ns(build_start);
                }
            }
        }));
    }
    for (auto& future : futures) {
        future.get();
    }
    std::vector<ClusterSubgraph> subgraphs;
    subgraphs.reserve(clusters.size());
    for (auto& slot : slots) {
        subgraphs.push_back(std::move(*slot));
    }
    return subgraphs;
}


struct ParentlessAdvanceAccounting {
    size_t raw_event_count = 0;
    size_t events_at_or_before_cutoff = 0;
    size_t events_after_cutoff = 0;
    cumulative_time_int stop_algorithmic_time = 0;
    // Source detector nodes that participated in MWPM notifications while this
    // cluster was advanced at the current level.  These are the active sources
    // touched by flooder events.  They can include active detectors imported
    // from lower levels, so this is intentionally not restricted to the
    // cluster's own local_active_detectors.
    std::vector<uint64_t> touched_active_detector_sources;
};

size_t parentless_env_size_t_or_zero(const char* name) {
    const char* value = std::getenv(name);
    if (value == nullptr || *value == 0) {
        return 0;
    }
    try {
        return static_cast<size_t>(std::stoull(value));
    } catch (...) {
        return 0;
    }
}

size_t parentless_node_id_from_detector_pointer(const Mwpm& mwpm, const DetectorNode* node) {
    if (node == nullptr || mwpm.flooder.graph.nodes.empty()) {
        return SIZE_MAX;
    }
    const auto* first = mwpm.flooder.graph.nodes.data();
    const auto* last = first + mwpm.flooder.graph.nodes.size();
    if (node < first || node >= last) {
        return SIZE_MAX;
    }
    return static_cast<size_t>(node - first);
}

void parentless_note_touched_source(
    const Mwpm& mwpm,
    const DetectorNode* node,
    std::vector<uint64_t>& touched_active_detector_sources) {
    const auto node_id = parentless_node_id_from_detector_pointer(mwpm, node);
    if (node_id == SIZE_MAX) {
        return;
    }
    const uint64_t source = static_cast<uint64_t>(node_id);
    if (std::find(
            touched_active_detector_sources.begin(),
            touched_active_detector_sources.end(),
            source) == touched_active_detector_sources.end()) {
        touched_active_detector_sources.push_back(source);
    }
}

void parentless_note_touched_sources_from_mwpm_event(
    const Mwpm& mwpm,
    const MwpmEvent& event,
    std::vector<uint64_t>& touched_active_detector_sources) {
    switch (event.event_type) {
        case REGION_HIT_REGION:
            parentless_note_touched_source(
                mwpm, event.region_hit_region_event_data.edge.loc_from, touched_active_detector_sources);
            parentless_note_touched_source(
                mwpm, event.region_hit_region_event_data.edge.loc_to, touched_active_detector_sources);
            break;
        case REGION_HIT_BOUNDARY:
            parentless_note_touched_source(
                mwpm, event.region_hit_boundary_event_data.edge.loc_from, touched_active_detector_sources);
            break;
        case BLOSSOM_SHATTER:
            // A shatter notification does not correspond to a new flooder touch
            // between active sources.  The active sources involved in creating
            // the blossom were already recorded by earlier hit-region events.
            break;
        case NO_EVENT:
            break;
        default:
            throw std::invalid_argument("Unrecognized MWPM event type while recording touched sources.");
    }
}

ParentlessAdvanceAccounting advance_initialized_mwpm_until_completion_counting_after(
    Mwpm& mwpm,
    cumulative_time_int cutoff_algorithmic_time,
    const std::string& debug_context = "") {
    ParentlessAdvanceAccounting out;
    const size_t mwpm_event_guard = parentless_env_size_t_or_zero("PYMATCHING_DEBUG_MAX_MWPM_EVENTS_PER_CLUSTER");
    const size_t progress_period = parentless_env_size_t_or_zero("PYMATCHING_DEBUG_PROGRESS_EVERY_MWPM_EVENTS");
    const std::string previous_debug_context = graph_flooder_debug_context;
    graph_flooder_debug_context = debug_context;
    try {
        while (true) {
            auto event = mwpm.flooder.run_until_next_mwpm_notification();
            if (event.event_type == pm::NO_EVENT) {
                break;
            }
            out.raw_event_count++;
            if (mwpm_event_guard != 0 && out.raw_event_count > mwpm_event_guard) {
                std::cerr
                    << "PYMATCHING DEBUG GUARD: exceeded "
                    << mwpm_event_guard
                    << " MWPM notifications in one parentless cluster"
                    << " context=" << debug_context
                    << " cur_time=" << mwpm.flooder.queue.cur_time
                    << " events_after_cutoff=" << out.events_after_cutoff
                    << " cutoff=" << cutoff_algorithmic_time
                    << " node_allocated=" << mwpm.node_arena.allocated.size()
                    << " node_available=" << mwpm.node_arena.available.size()
                    << std::endl;
                throw std::runtime_error(
                    "PYMATCHING_DEBUG_MAX_MWPM_EVENTS_PER_CLUSTER exceeded; see stderr for context");
            }
            if (progress_period != 0 && out.raw_event_count % progress_period == 0) {
                std::cerr
                    << "PYMATCHING DEBUG PROGRESS:"
                    << " context=" << debug_context
                    << " raw_events=" << out.raw_event_count
                    << " events_after_cutoff=" << out.events_after_cutoff
                    << " cur_time=" << mwpm.flooder.queue.cur_time
                    << " cutoff=" << cutoff_algorithmic_time
                    << std::endl;
            }
            // run_until_next_mwpm_notification advances the flooder queue clock to
            // the algorithmic time of the returned MWPM notification. Events at or
            // before the previous level's maximum stop time are considered to have
            // overlapped lower-level cluster execution, so only strictly later
            // notifications extend the level-by-level critical path.
            if (mwpm.flooder.queue.cur_time > cutoff_algorithmic_time) {
                out.events_after_cutoff++;
            } else {
                out.events_at_or_before_cutoff++;
            }
            parentless_note_touched_sources_from_mwpm_event(
                mwpm, event, out.touched_active_detector_sources);
            mwpm.process_event(event);
        }
        if (mwpm.node_arena.allocated.size() != mwpm.node_arena.available.size()) {
            mwpm.reset();
            throw std::invalid_argument(
                "No perfect matching could be found while advancing an initialized parentless pipeline cluster.");
        }
        out.stop_algorithmic_time = mwpm.flooder.queue.cur_time;
    } catch (...) {
        graph_flooder_debug_context = previous_debug_context;
        throw;
    }
    graph_flooder_debug_context = previous_debug_context;
    return out;
}

size_t advance_initialized_mwpm_until_completion(Mwpm& mwpm) {
    return advance_initialized_mwpm_until_completion_counting_after(
        mwpm,
        std::numeric_limits<cumulative_time_int>::min(),
        "legacy_completion").raw_event_count;
}

void parentless_pipeline_release_cluster_mwpm_if_needed(ClusterSubgraph& subgraph) {
    if (subgraph.persistent_mwpm != nullptr) {
        release_persistent_full_graph_worker_mwpm(subgraph.persistent_mwpm);
        subgraph.persistent_mwpm = nullptr;
        subgraph.persistent_mwpm_acquired_clean = false;
    }
}

void parentless_pipeline_ensure_cluster_initialized(
    UserGraph& graph,
    LockstepDecodeResult& result,
    size_t cluster_id,
    const LockstepSchedulerConfig& config,
    std::vector<uint8_t>& cluster_mwpm_initialized) {
    if (cluster_id < cluster_mwpm_initialized.size() && cluster_mwpm_initialized[cluster_id]) {
        return;
    }
    auto cluster_init_start = steady_clock::now();
    auto& subgraph = result.cluster_subgraphs[cluster_id];
    result.cluster_extraction_detection_events[cluster_id] = subgraph.local_active_detectors;
    rebuild_extraction_seen_for_cluster(result, cluster_id);
    acquire_persistent_full_graph_worker_mwpm_for_subgraph(graph, subgraph);
    auto& mwpm = cluster_mwpm(subgraph, config.edge_correlations);
    if (subgraph.persistent_mwpm != nullptr && !subgraph.persistent_mwpm_acquired_clean) {
        mwpm.reset();
    }
    initialize_mwpm_for_detection_events(mwpm, subgraph.local_active_detectors);
    subgraph.persistent_mwpm_acquired_clean = false;
    if (cluster_id < result.cluster_initial_setup_wall_ns.size()) {
        result.cluster_initial_setup_wall_ns[cluster_id] += elapsed_ns(cluster_init_start);
    }
    if (cluster_id < cluster_mwpm_initialized.size()) {
        cluster_mwpm_initialized[cluster_id] = 1;
    }
}

std::vector<uint64_t> parentless_pipeline_observable_extraction_sources(
    LockstepDecodeResult& result,
    size_t cluster_id,
    const std::vector<size_t>& imported_source_cluster_ids,
    const std::vector<uint64_t>& touched_active_detector_sources,
    std::vector<uint8_t>& direct_result_superseded) {
    if (cluster_id >= result.cluster_subgraphs.size()) {
        throw std::invalid_argument("Parentless observable extraction cluster id is outside the subgraph table.");
    }

    const auto& target_subgraph = result.cluster_subgraphs[cluster_id];
    std::vector<uint64_t> extraction_sources = target_subgraph.local_active_detectors;

    auto source_was_touched = [&](const std::vector<uint64_t>& remapped_sources) {
        for (auto source : remapped_sources) {
            if (std::find(
                    touched_active_detector_sources.begin(),
                    touched_active_detector_sources.end(),
                    source) != touched_active_detector_sources.end()) {
                return true;
            }
        }
        return false;
    };

    // An untouched inherited matching configuration keeps its lower-level
    // direct result as its terminal contribution.  If the current cluster
    // touches any detector belonging to an inherited configuration, the whole
    // inherited configuration is replaced by the current cluster's updated
    // configuration.  Include all detector sources of that configuration in
    // this cluster's observable extraction and suppress the stale lower-level
    // direct result during final aggregation.
    for (auto source_cluster_id : imported_source_cluster_ids) {
        if (source_cluster_id >= result.cluster_subgraphs.size() ||
            source_cluster_id >= result.cluster_stopped_states.size()) {
            throw std::invalid_argument("Parentless observable extraction source cluster id is outside saved state tables.");
        }
        auto remapped_sources = remap_source_detection_events_from_stopped_state(
            target_subgraph,
            result.cluster_subgraphs[source_cluster_id],
            result.cluster_stopped_states[source_cluster_id]);
        if (!source_was_touched(remapped_sources)) {
            continue;
        }
        if (source_cluster_id >= direct_result_superseded.size()) {
            direct_result_superseded.resize(result.hierarchy.clusters.size(), 0);
        }
        direct_result_superseded[source_cluster_id] = 1;
        extraction_sources.insert(
            extraction_sources.end(), remapped_sources.begin(), remapped_sources.end());
    }

    // Keep the event-level touched sources themselves as a safety net.  Own
    // active detectors are already present, while imported touched detectors
    // may have come from a configuration whose stopped-state representation is
    // intentionally sparse.
    extraction_sources.insert(
        extraction_sources.end(),
        touched_active_detector_sources.begin(),
        touched_active_detector_sources.end());
    std::sort(extraction_sources.begin(), extraction_sources.end());
    extraction_sources.erase(
        std::unique(extraction_sources.begin(), extraction_sources.end()),
        extraction_sources.end());
    return extraction_sources;
}

void parentless_pipeline_export_or_extract_stopped_cluster(
    LockstepDecodeResult& result,
    size_t cluster_id,
    const LockstepSchedulerConfig& config,
    const std::vector<size_t>& imported_source_cluster_ids,
    const std::vector<uint64_t>& touched_active_detector_sources,
    std::vector<uint8_t>& direct_result_superseded,
    bool force_export_stopped_state = false) {
    auto& subgraph = result.cluster_subgraphs[cluster_id];
    auto& mwpm = cluster_mwpm(subgraph, config.edge_correlations);
    const bool is_root_cluster =
        result.hierarchy.clusters[cluster_id].parent_id == NO_PROCESSING_CLUSTER_PARENT;

    result.cluster_execution_states[cluster_id] = LockstepClusterExecutionState::QUIESCENT_STOPPED;
    result.cluster_has_stopped_state[cluster_id] = true;

    const auto observable_extraction_sources = parentless_pipeline_observable_extraction_sources(
        result,
        cluster_id,
        imported_source_cluster_ids,
        touched_active_detector_sources,
        direct_result_superseded);

    // In the parentless level pipeline, root/non-root labels are only an
    // aggregation detail.  A cluster that has a next-level broadcast target must
    // export its stopped state even if the hierarchy builder marked it as a root
    // (this happens when the residual syndrome is exhausted and no empty terminal
    // root is appended).  Otherwise lower-level states would not be inherited by
    // level+1 clusters and the critical-event pipeline would no longer match the
    // intended broadcast model.
    if (is_root_cluster && !force_export_stopped_state) {
        const bool already_has_direct_result =
            (cluster_id < result.cluster_has_direct_compact_result.size() &&
             result.cluster_has_direct_compact_result[cluster_id]) ||
            (cluster_id < result.cluster_has_direct_result.size() && result.cluster_has_direct_result[cluster_id]);
        if (already_has_direct_result) {
            result.cluster_result_is_provisional[cluster_id] = false;
            return;
        }
        auto extract_start = steady_clock::now();
        if (mwpm.flooder.graph.num_observables <= sizeof(pm::obs_int) * 8) {
            result.cluster_direct_compact_results[cluster_id] = extract_compact_result_from_current_mwpm_state(
                mwpm, observable_extraction_sources);
            result.cluster_has_direct_compact_result[cluster_id] = true;
            result.cluster_debug_compact_results[cluster_id] = result.cluster_direct_compact_results[cluster_id];
            result.cluster_has_debug_compact_result[cluster_id] = true;
        } else {
            result.cluster_direct_results[cluster_id] = extract_result_from_current_mwpm_state(
                mwpm, observable_extraction_sources);
            result.cluster_has_direct_result[cluster_id] = true;
        }
        const auto elapsed = elapsed_ns(extract_start);
        result.profiling_stats.root_extraction_wall_ns += elapsed;
        if (cluster_id < result.cluster_mark_stop_wall_ns.size()) {
            result.cluster_mark_stop_wall_ns[cluster_id] += elapsed;
        }
        if (subgraph.persistent_mwpm != nullptr) {
            parentless_pipeline_release_cluster_mwpm_if_needed(subgraph);
        } else if (subgraph.use_shared_full_graph_mwpm) {
            mwpm.reset();
        }
        result.cluster_result_is_provisional[cluster_id] = false;
        return;
    }

    auto export_start = steady_clock::now();
    result.cluster_stopped_states[cluster_id] = export_mwpm_stopped_state_for_active_detectors(
        mwpm, touched_active_detector_sources);
    if (mwpm.flooder.graph.num_observables <= sizeof(pm::obs_int) * 8) {
        result.cluster_direct_compact_results[cluster_id] = extract_compact_result_from_current_mwpm_state(
            mwpm, observable_extraction_sources);
        result.cluster_has_direct_compact_result[cluster_id] = true;
        result.cluster_debug_compact_results[cluster_id] = result.cluster_direct_compact_results[cluster_id];
        result.cluster_has_debug_compact_result[cluster_id] = true;
    } else {
        result.cluster_direct_results[cluster_id] = extract_result_from_current_mwpm_state(
            mwpm, observable_extraction_sources);
        result.cluster_has_direct_result[cluster_id] = true;
    }
    const auto elapsed = elapsed_ns(export_start);
    if (cluster_id < result.cluster_mark_stop_wall_ns.size()) {
        result.cluster_mark_stop_wall_ns[cluster_id] += elapsed;
    }
    parentless_pipeline_release_cluster_mwpm_if_needed(subgraph);
    result.cluster_result_is_provisional[cluster_id] = true;
}

std::vector<size_t> parentless_all_upper_targets_for_child(
    const LockstepDecodeResult& result,
    size_t child_cluster_id) {
    if (child_cluster_id >= result.hierarchy.clusters.size()) {
        return {};
    }
    const auto child_level = result.hierarchy.clusters[child_cluster_id].level;
    std::vector<size_t> targets;
    for (size_t target_level = child_level + 1;
         target_level < result.hierarchy.cluster_ids_by_level.size();
         target_level++) {
        for (auto target_cluster_id : result.hierarchy.cluster_ids_by_level[target_level]) {
            if (target_cluster_id != child_cluster_id) {
                targets.push_back(target_cluster_id);
            }
        }
    }
    return targets;
}

bool parentless_pipeline_child_has_upper_targets(
    const LockstepDecodeResult& result,
    size_t child_cluster_id) {
    return !parentless_all_upper_targets_for_child(result, child_cluster_id).empty();
}

std::vector<size_t> parentless_empty_terminal_roots(const LockstepDecodeResult& result) {
    std::vector<size_t> roots;
    for (auto root_cluster_id : result.hierarchy.root_cluster_ids) {
        if (cluster_is_empty_terminal_root(result, root_cluster_id)) {
            roots.push_back(root_cluster_id);
        }
    }
    return roots;
}

void parentless_pipeline_import_sources_into_target(
    UserGraph& graph,
    LockstepDecodeResult& result,
    size_t target_cluster_id,
    const std::vector<size_t>& source_cluster_ids,
    const LockstepSchedulerConfig& config,
    std::vector<uint8_t>& cluster_mwpm_initialized,
    std::vector<std::vector<size_t>>& imported_source_cluster_ids_by_target) {
    if (source_cluster_ids.empty()) {
        return;
    }
    auto import_start = steady_clock::now();
    parentless_pipeline_ensure_cluster_initialized(
        graph, result, target_cluster_id, config, cluster_mwpm_initialized);
    if (cluster_is_empty_terminal_root(result, target_cluster_id)) {
        terminal_import_children_into_empty_root(result, target_cluster_id, source_cluster_ids, config);
    } else {
        auto& target_subgraph = result.cluster_subgraphs[target_cluster_id];
        std::vector<size_t> imported_target_local_node_ids;
        for (auto source_cluster_id : source_cluster_ids) {
            StoppedStateImportOptions options;
            options.reschedule_after_import = false;
            options.imported_target_local_node_ids = &imported_target_local_node_ids;
            options.skip_unmappable_or_conflicting_state = false;
            options.allow_detector_ownership_overwrite = true;
            import_stopped_state_into_cluster_subgraph(
                target_subgraph,
                result.cluster_subgraphs[source_cluster_id],
                result.cluster_stopped_states[source_cluster_id],
                options,
                config.edge_correlations);
            result.profiling_stats.child_import_count++;
        }
        reschedule_imported_detector_nodes(target_subgraph, imported_target_local_node_ids, config.edge_correlations);
        if (!imported_target_local_node_ids.empty() &&
            result.cluster_execution_states[target_cluster_id] != LockstepClusterExecutionState::RUNNING) {
            result.cluster_execution_states[target_cluster_id] = LockstepClusterExecutionState::RUNNING;
            result.cluster_has_stopped_state[target_cluster_id] = false;
            result.cluster_result_is_provisional[target_cluster_id] = false;
            if (target_cluster_id < result.cluster_has_direct_compact_result.size()) {
                result.cluster_has_direct_compact_result[target_cluster_id] = false;
            }
            if (target_cluster_id < result.cluster_has_direct_result.size()) {
                result.cluster_has_direct_result[target_cluster_id] = false;
            }
        }
    }
    if (target_cluster_id >= imported_source_cluster_ids_by_target.size()) {
        throw std::invalid_argument("Parentless import target is outside imported-source tracking table.");
    }
    auto& imported_sources = imported_source_cluster_ids_by_target[target_cluster_id];
    for (auto source_cluster_id : source_cluster_ids) {
        if (std::find(imported_sources.begin(), imported_sources.end(), source_cluster_id) == imported_sources.end()) {
            imported_sources.push_back(source_cluster_id);
        }
    }

    result.profiling_stats.scheduler_batch_count++;
    result.profiling_stats.import_candidate_count += source_cluster_ids.size();
    const auto import_elapsed = elapsed_ns(import_start);
    result.profiling_stats.child_import_wall_ns += import_elapsed;
    if (target_cluster_id < result.cluster_child_import_wall_ns.size()) {
        result.cluster_child_import_wall_ns[target_cluster_id] += import_elapsed;
    }
}

void parentless_pipeline_apply_pending_imports_for_level(
    UserGraph& graph,
    LockstepDecodeResult& result,
    size_t target_level,
    std::vector<std::vector<size_t>>& pending_import_sources_by_target,
    const LockstepSchedulerConfig& config,
    std::vector<uint8_t>& cluster_mwpm_initialized,
    std::vector<std::vector<size_t>>& imported_source_cluster_ids_by_target) {
    if (target_level >= result.hierarchy.cluster_ids_by_level.size()) {
        return;
    }
    for (auto target_cluster_id : result.hierarchy.cluster_ids_by_level[target_level]) {
        if (target_cluster_id >= pending_import_sources_by_target.size()) {
            continue;
        }
        auto& pending = pending_import_sources_by_target[target_cluster_id];
        if (pending.empty()) {
            continue;
        }
        parentless_pipeline_import_sources_into_target(
            graph,
            result,
            target_cluster_id,
            pending,
            config,
            cluster_mwpm_initialized,
            imported_source_cluster_ids_by_target);
        pending.clear();
    }
}

void parentless_pipeline_queue_finished_level_imports(
    LockstepDecodeResult& result,
    size_t child_level,
    std::vector<std::vector<size_t>>& pending_import_sources_by_target) {
    if (child_level >= result.hierarchy.cluster_ids_by_level.size()) {
        return;
    }
    std::vector<size_t> queued_source_cluster_ids;
    for (auto source_cluster_id : result.hierarchy.cluster_ids_by_level[child_level]) {
        if (!result.cluster_has_stopped_state[source_cluster_id]) {
            continue;
        }
        if (result.cluster_execution_states[source_cluster_id] == LockstepClusterExecutionState::IMPORTED_INTO_PARENT) {
            continue;
        }
        const auto targets = parentless_all_upper_targets_for_child(result, source_cluster_id);
        if (targets.empty()) {
            continue;
        }
        for (auto target_cluster_id : targets) {
            if (target_cluster_id >= pending_import_sources_by_target.size()) {
                continue;
            }
            auto& pending = pending_import_sources_by_target[target_cluster_id];
            if (std::find(pending.begin(), pending.end(), source_cluster_id) == pending.end()) {
                pending.push_back(source_cluster_id);
            }
        }
        queued_source_cluster_ids.push_back(source_cluster_id);
    }
    for (auto source_cluster_id : queued_source_cluster_ids) {
        result.cluster_execution_states[source_cluster_id] = LockstepClusterExecutionState::IMPORTED_INTO_PARENT;
        result.cluster_result_is_provisional[source_cluster_id] = false;
    }
}

void aggregate_parentless_pipeline_roots(
    UserGraph& graph,
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config,
    const std::vector<uint8_t>& direct_result_superseded) {
    const bool can_use_compact_root_aggregation = graph.get_num_observables() <= sizeof(pm::obs_int) * 8;
    pm::MatchingResult compact_root_aggregate;
    std::vector<uint8_t> seen_root_ids(result.hierarchy.clusters.size(), 0);
    for (auto root_cluster_id : result.hierarchy.root_cluster_ids) {
        if (root_cluster_id >= result.hierarchy.clusters.size()) {
            throw std::invalid_argument("Parentless pipeline root id is outside the cluster hierarchy.");
        }
        if (seen_root_ids[root_cluster_id]) {
            continue;
        }
        seen_root_ids[root_cluster_id] = 1;
        if (root_cluster_id < direct_result_superseded.size() && direct_result_superseded[root_cluster_id]) {
            continue;
        }
        const bool root_imported_into_upper_cluster =
            root_cluster_id < result.cluster_execution_states.size() &&
            result.cluster_execution_states[root_cluster_id] == LockstepClusterExecutionState::IMPORTED_INTO_PARENT;
        const bool root_has_saved_own_result =
            (root_cluster_id < result.cluster_has_direct_compact_result.size() &&
             result.cluster_has_direct_compact_result[root_cluster_id]) ||
            (root_cluster_id < result.cluster_has_direct_result.size() &&
             result.cluster_has_direct_result[root_cluster_id]);
        if (root_imported_into_upper_cluster && !root_has_saved_own_result) {
            continue;
        }
        if (can_use_compact_root_aggregation) {
            if (root_cluster_id < result.cluster_has_direct_compact_result.size() &&
                result.cluster_has_direct_compact_result[root_cluster_id]) {
                compact_root_aggregate += result.cluster_direct_compact_results[root_cluster_id];
                continue;
            }
            auto& root_subgraph = result.cluster_subgraphs[root_cluster_id];
            auto& root_mwpm = cluster_mwpm(root_subgraph, config.edge_correlations);
            auto extract_start = steady_clock::now();
            compact_root_aggregate += extract_compact_result_from_current_mwpm_state(
                root_mwpm, result.cluster_extraction_detection_events[root_cluster_id]);
            result.profiling_stats.root_extraction_wall_ns += elapsed_ns(extract_start);
            if (root_subgraph.persistent_mwpm != nullptr) {
                release_persistent_full_graph_worker_mwpm(root_subgraph.persistent_mwpm);
                root_subgraph.persistent_mwpm = nullptr;
            } else if (root_subgraph.use_shared_full_graph_mwpm) {
                root_mwpm.reset();
            }
            continue;
        }
        if (root_cluster_id < result.cluster_has_direct_result.size() &&
            result.cluster_has_direct_result[root_cluster_id]) {
            result.root_aggregate_result += result.cluster_direct_results[root_cluster_id];
            continue;
        }
        auto& root_subgraph = result.cluster_subgraphs[root_cluster_id];
        auto& root_mwpm = cluster_mwpm(root_subgraph, config.edge_correlations);
        auto extract_start = steady_clock::now();
        result.root_aggregate_result += extract_result_from_current_mwpm_state(
            root_mwpm, result.cluster_extraction_detection_events[root_cluster_id]);
        result.profiling_stats.root_extraction_wall_ns += elapsed_ns(extract_start);
        if (root_subgraph.persistent_mwpm != nullptr) {
            release_persistent_full_graph_worker_mwpm(root_subgraph.persistent_mwpm);
            root_subgraph.persistent_mwpm = nullptr;
        } else if (root_subgraph.use_shared_full_graph_mwpm) {
            root_mwpm.reset();
        }
    }
    if (can_use_compact_root_aggregation) {
        fill_bit_vector_from_obs_mask(
            compact_root_aggregate.obs_mask,
            result.root_aggregate_result.obs_crossed.data(),
            graph.get_num_observables());
        result.root_aggregate_result.weight = compact_root_aggregate.weight;
    }
}

bool parentless_pipeline_all_clusters_are_roots(const LockstepDecodeResult& result) {
    if (result.hierarchy.clusters.empty()) {
        return false;
    }
    for (const auto& cluster : result.hierarchy.clusters) {
        if (cluster.parent_id != NO_PROCESSING_CLUSTER_PARENT) {
            return false;
        }
    }
    return result.hierarchy.root_cluster_ids.size() == result.hierarchy.clusters.size();
}

bool run_parentless_all_root_direct_io_path(
    UserGraph& graph,
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config) {
    if (!config.enable_full_graph_worker_subgraphs) {
        return false;
    }
    if (config.edge_correlations) {
        return false;
    }
    if (graph.get_num_observables() > sizeof(pm::obs_int) * 8) {
        return false;
    }
    if (!config.enable_parentless_all_root_direct_io_path) {
        return false;
    }
    if (!parentless_pipeline_all_clusters_are_roots(result)) {
        return false;
    }

    const auto& root_ids = result.hierarchy.root_cluster_ids;
    std::vector<pm::MatchingResult> root_results(root_ids.size());

    auto scheduler_start = steady_clock::now();
    const bool enable_actual_parallel_level_execution =
        std::getenv("PYMATCHING_PARENTLESS_ACTUAL_PARALLEL_LEVEL") != nullptr &&
        std::string(std::getenv("PYMATCHING_PARENTLESS_ACTUAL_PARALLEL_LEVEL")) == "1";
    const size_t level_workers = std::min(checked_num_workers(config), root_ids.size());

    auto run_one_root = [&](size_t k) {
        const size_t cluster_id = root_ids[k];
        auto& subgraph = result.cluster_subgraphs[cluster_id];
        acquire_persistent_full_graph_worker_mwpm_for_subgraph(graph, subgraph);
        auto& mwpm = cluster_mwpm(subgraph, config.edge_correlations);
        auto step_start = steady_clock::now();
        uint64_t init_elapsed = 0;
        uint64_t advance_elapsed = 0;
        uint64_t extract_elapsed = 0;
        uint64_t release_elapsed = 0;
        uint64_t event_count = 0;
        try {
            // Direct all-root I/O split into explicit stages. This measures the
            // same worker operation as decode_detection_events_for_up_to_64_observables(...),
            // but exposes whether the cost is active-detector insertion, sparse-blossom
            // timeline advance, compact observable extraction, or state release/reset.
            auto init_start = steady_clock::now();
            initialize_mwpm_for_detection_events(mwpm, subgraph.local_active_detectors);
            init_elapsed = elapsed_ns(init_start);

            auto advance_start = steady_clock::now();
            while (true) {
                auto event = mwpm.flooder.run_until_next_mwpm_notification();
                if (event.event_type == pm::NO_EVENT) {
                    break;
                }
                event_count++;
                mwpm.process_event(event);
            }
            if (mwpm.node_arena.allocated.size() != mwpm.node_arena.available.size()) {
                mwpm.reset();
                throw std::invalid_argument(
                    "No perfect matching could be found in direct all-root worker advance.");
            }
            advance_elapsed = elapsed_ns(advance_start);
            const auto stop_algorithmic_time = mwpm.flooder.queue.cur_time;
            if (cluster_id < result.cluster_critical_path_event_counts.size()) {
                result.cluster_critical_path_event_counts[cluster_id] = event_count;
            }
            if (cluster_id < result.cluster_has_critical_path_event_count.size()) {
                result.cluster_has_critical_path_event_count[cluster_id] = 1;
            }
            if (cluster_id < result.cluster_stop_algorithmic_times.size()) {
                result.cluster_stop_algorithmic_times[cluster_id] = stop_algorithmic_time;
            }
            if (cluster_id < result.cluster_critical_path_algorithmic_times.size()) {
                result.cluster_critical_path_algorithmic_times[cluster_id] = stop_algorithmic_time;
            }
            if (cluster_id < result.cluster_has_cutoff_algorithmic_time.size()) {
                result.cluster_has_cutoff_algorithmic_time[cluster_id] = 0;
            }
            if (cluster_id < result.cluster_cutoff_algorithmic_times.size()) {
                result.cluster_cutoff_algorithmic_times[cluster_id] = 0;
            }
            if (cluster_id < result.cluster_events_at_or_before_cutoff.size()) {
                result.cluster_events_at_or_before_cutoff[cluster_id] = 0;
            }
            if (cluster_id < result.cluster_events_after_cutoff.size()) {
                result.cluster_events_after_cutoff[cluster_id] = event_count;
            }

            auto extract_start = steady_clock::now();
            root_results[k] = extract_compact_result_from_current_mwpm_state(
                mwpm, subgraph.local_active_detectors);
            extract_elapsed = elapsed_ns(extract_start);
        } catch (...) {
            auto release_start = steady_clock::now();
            if (subgraph.persistent_mwpm != nullptr) {
                parentless_pipeline_release_cluster_mwpm_if_needed(subgraph);
            } else if (subgraph.use_shared_full_graph_mwpm) {
                mwpm.reset();
            }
            release_elapsed = elapsed_ns(release_start);
            throw;
        }
        auto release_start = steady_clock::now();
        if (subgraph.persistent_mwpm != nullptr) {
            parentless_pipeline_release_cluster_mwpm_if_needed(subgraph);
        } else if (subgraph.use_shared_full_graph_mwpm) {
            mwpm.reset();
        }
        release_elapsed = elapsed_ns(release_start);
        const auto step_elapsed = elapsed_ns(step_start);
        if (cluster_id < result.cluster_step_wall_ns.size()) {
            result.cluster_step_wall_ns[cluster_id] += step_elapsed;
        }
        if (cluster_id < result.cluster_step_event_counts.size()) {
            result.cluster_step_event_counts[cluster_id] += event_count;
        }
        auto& stats = result.profiling_stats;
        stats.direct_worker_init_wall_ns += init_elapsed;
        stats.direct_worker_advance_wall_ns += advance_elapsed;
        stats.direct_worker_extract_wall_ns += extract_elapsed;
        stats.direct_worker_release_wall_ns += release_elapsed;
        stats.direct_worker_total_wall_ns += step_elapsed;
        stats.max_direct_worker_init_wall_ns = std::max(stats.max_direct_worker_init_wall_ns, init_elapsed);
        stats.max_direct_worker_advance_wall_ns = std::max(stats.max_direct_worker_advance_wall_ns, advance_elapsed);
        stats.max_direct_worker_extract_wall_ns = std::max(stats.max_direct_worker_extract_wall_ns, extract_elapsed);
        stats.max_direct_worker_release_wall_ns = std::max(stats.max_direct_worker_release_wall_ns, release_elapsed);
        stats.max_direct_worker_total_wall_ns = std::max(stats.max_direct_worker_total_wall_ns, step_elapsed);
        stats.direct_worker_active_detectors += subgraph.local_active_detectors.size();
        stats.max_direct_worker_active_detectors = std::max<uint64_t>(
            stats.max_direct_worker_active_detectors, subgraph.local_active_detectors.size());
        result.cluster_execution_states[cluster_id] = LockstepClusterExecutionState::QUIESCENT_STOPPED;
        result.cluster_has_stopped_state[cluster_id] = true;
        result.cluster_has_direct_compact_result[cluster_id] = true;
        result.cluster_direct_compact_results[cluster_id] = root_results[k];
        result.cluster_result_is_provisional[cluster_id] = false;
    };

    if (enable_actual_parallel_level_execution && level_workers > 1) {
        auto& pool = parentless_level_worker_pool(level_workers);
        pool.parallel_for(root_ids.size(), run_one_root);
    } else {
        for (size_t k = 0; k < root_ids.size(); k++) {
            run_one_root(k);
        }
    }
    result.profiling_stats.scheduler_cluster_step_count += root_ids.size();
    result.profiling_stats.scheduler_wall_ns += elapsed_ns(scheduler_start);

    pm::MatchingResult compact_root_aggregate;
    for (const auto& r : root_results) {
        compact_root_aggregate += r;
    }
    fill_bit_vector_from_obs_mask(
        compact_root_aggregate.obs_mask,
        result.root_aggregate_result.obs_crossed.data(),
        graph.get_num_observables());
    result.root_aggregate_result.weight = compact_root_aggregate.weight;

    finalize_cluster_step_parallel_accounting(result, config);
    return true;
}

bool run_parentless_level_pipeline(
    UserGraph& graph,
    LockstepDecodeResult& result,
    const LockstepSchedulerConfig& config) {
    if (!config.enable_parentless_level_pipeline) {
        return false;
    }
    if (!config.enable_parentless_broadcast_import || !config.enable_growing_only_processing_clusters) {
        return false;
    }
    if (config.edge_correlations) {
        throw std::invalid_argument("Parentless level pipeline does not yet support edge correlations.");
    }

    if (run_parentless_all_root_direct_io_path(graph, result, config)) {
        return true;
    }

    // In the all-upper-level parentless broadcast model, every processing
    // cluster owns a disjoint set of original active detectors.  Root labels from
    // the hierarchy builder are not a parent/child aggregation contract here;
    // aggregate each cluster's own active-detector contribution exactly once.
    result.hierarchy.root_cluster_ids.clear();
    result.hierarchy.root_cluster_ids.reserve(result.hierarchy.clusters.size());
    for (const auto& cluster : result.hierarchy.clusters) {
        result.hierarchy.root_cluster_ids.push_back(cluster.id);
    }
    result.profiling_stats.root_cluster_count = result.hierarchy.root_cluster_ids.size();

    auto setup_start = steady_clock::now();
    for (size_t cluster_id = 0; cluster_id < result.cluster_subgraphs.size(); cluster_id++) {
        result.cluster_extraction_detection_events[cluster_id] =
            result.cluster_subgraphs[cluster_id].local_active_detectors;
        rebuild_extraction_seen_for_cluster(result, cluster_id);
    }
    std::vector<uint8_t> cluster_mwpm_initialized(result.cluster_subgraphs.size(), 0);
    result.profiling_stats.initial_mwpm_setup_wall_ns += elapsed_ns(setup_start);

    auto scheduler_start = steady_clock::now();
    size_t prior_level_critical_events = 0;
    cumulative_time_int prior_level_stop_time = 0;
    bool have_prior_completed_level = false;
    std::vector<std::vector<size_t>> pending_import_sources_by_target(result.hierarchy.clusters.size());
    std::vector<std::vector<size_t>> imported_source_cluster_ids_by_target(result.hierarchy.clusters.size());
    std::vector<uint8_t> direct_result_superseded(result.hierarchy.clusters.size(), 0);
    for (size_t level = 0; level < result.hierarchy.cluster_ids_by_level.size(); level++) {
        parentless_pipeline_apply_pending_imports_for_level(
            graph,
            result,
            level,
            pending_import_sources_by_target,
            config,
            cluster_mwpm_initialized,
            imported_source_cluster_ids_by_target);
        const auto& cluster_ids = result.hierarchy.cluster_ids_by_level[level];
        std::vector<size_t> runnable_cluster_ids;
        runnable_cluster_ids.reserve(cluster_ids.size());
        for (auto cluster_id : cluster_ids) {
            if (result.cluster_execution_states[cluster_id] == LockstepClusterExecutionState::IMPORTED_INTO_PARENT) {
                continue;
            }
            const bool already_has_root_result =
                (cluster_id < result.cluster_has_direct_compact_result.size() &&
                 result.cluster_has_direct_compact_result[cluster_id]) ||
                (cluster_id < result.cluster_has_direct_result.size() && result.cluster_has_direct_result[cluster_id]);
            if (already_has_root_result &&
                result.hierarchy.clusters[cluster_id].parent_id == NO_PROCESSING_CLUSTER_PARENT) {
                continue;
            }
            runnable_cluster_ids.push_back(cluster_id);
        }

        // Safety-first accounting path: actually execute the sparse-blossom worker
        // state for every cluster, but do it sequentially. The event metric models
        // ideal level overlap by allowing current-level notifications at or before
        // the previous level's maximum stop time to run in parallel with the lower
        // levels.  The overlap phase costs max(lower-level critical events,
        // current-cluster events at/before cutoff), and notifications strictly
        // after the cutoff are then added.
        const cumulative_time_int cutoff_time = have_prior_completed_level
            ? prior_level_stop_time
            : std::numeric_limits<cumulative_time_int>::min();
        size_t level_max_critical_events = prior_level_critical_events;
        cumulative_time_int level_max_stop_time = prior_level_stop_time;
        for (auto cluster_id : runnable_cluster_ids) {
            parentless_pipeline_ensure_cluster_initialized(
                graph, result, cluster_id, config, cluster_mwpm_initialized);
            auto& mwpm = cluster_mwpm(result.cluster_subgraphs[cluster_id], config.edge_correlations);
            auto step_start = steady_clock::now();
            const std::string debug_context =
                "parentless_level_pipeline"
                " level=" + std::to_string(level) +
                " cluster=" + std::to_string(cluster_id) +
                " parent=" + std::to_string(result.hierarchy.clusters[cluster_id].parent_id) +
                " active=" + std::to_string(result.cluster_subgraphs[cluster_id].local_active_detectors.size()) +
                " extraction_active=" + std::to_string(result.cluster_extraction_detection_events[cluster_id].size()) +
                " cutoff=" + std::to_string(cutoff_time);
            if (std::getenv("PYMATCHING_DEBUG_TRACE_PARENTLESS_LEVELS") != nullptr) {
                std::cerr << "PYMATCHING DEBUG START: " << debug_context << std::endl;
            }
            const auto advance = advance_initialized_mwpm_until_completion_counting_after(mwpm, cutoff_time, debug_context);
            if (std::getenv("PYMATCHING_DEBUG_TRACE_PARENTLESS_LEVELS") != nullptr) {
                std::cerr
                    << "PYMATCHING DEBUG DONE: " << debug_context
                    << " raw_events=" << advance.raw_event_count
                    << " events_after_cutoff=" << advance.events_after_cutoff
                    << " stop_time=" << advance.stop_algorithmic_time
                    << std::endl;
            }
            const auto step_elapsed = elapsed_ns(step_start);
            if (cluster_id < result.cluster_step_wall_ns.size()) {
                result.cluster_step_wall_ns[cluster_id] += step_elapsed;
            }
            if (cluster_id < result.cluster_step_event_counts.size()) {
                result.cluster_step_event_counts[cluster_id] += advance.raw_event_count;
            }
            const size_t cluster_overlap_event_base = std::max(
                prior_level_critical_events,
                advance.events_at_or_before_cutoff);
            const size_t cluster_critical_events =
                cluster_overlap_event_base + advance.events_after_cutoff;
            const cumulative_time_int cluster_critical_time = std::max(
                prior_level_stop_time,
                advance.stop_algorithmic_time);
            if (cluster_id < result.cluster_critical_path_event_counts.size()) {
                result.cluster_critical_path_event_counts[cluster_id] = cluster_critical_events;
            }
            if (cluster_id < result.cluster_has_critical_path_event_count.size()) {
                result.cluster_has_critical_path_event_count[cluster_id] = 1;
            }
            if (cluster_id < result.cluster_stop_algorithmic_times.size()) {
                result.cluster_stop_algorithmic_times[cluster_id] = advance.stop_algorithmic_time;
            }
            if (cluster_id < result.cluster_critical_path_algorithmic_times.size()) {
                result.cluster_critical_path_algorithmic_times[cluster_id] = cluster_critical_time;
            }
            if (cluster_id < result.cluster_has_cutoff_algorithmic_time.size()) {
                result.cluster_has_cutoff_algorithmic_time[cluster_id] = have_prior_completed_level ? 1 : 0;
            }
            if (cluster_id < result.cluster_cutoff_algorithmic_times.size()) {
                result.cluster_cutoff_algorithmic_times[cluster_id] = have_prior_completed_level ? cutoff_time : 0;
            }
            if (cluster_id < result.cluster_prior_level_critical_event_counts.size()) {
                result.cluster_prior_level_critical_event_counts[cluster_id] = prior_level_critical_events;
            }
            if (cluster_id < result.cluster_overlap_event_bases.size()) {
                result.cluster_overlap_event_bases[cluster_id] = cluster_overlap_event_base;
            }
            if (cluster_id < result.cluster_events_at_or_before_cutoff.size()) {
                result.cluster_events_at_or_before_cutoff[cluster_id] = advance.events_at_or_before_cutoff;
            }
            if (cluster_id < result.cluster_events_after_cutoff.size()) {
                result.cluster_events_after_cutoff[cluster_id] = advance.events_after_cutoff;
            }
            level_max_critical_events = std::max(level_max_critical_events, cluster_critical_events);
            level_max_stop_time = std::max(level_max_stop_time, advance.stop_algorithmic_time);
            result.profiling_stats.scheduler_cluster_step_count++;
            const bool has_upper_targets =
                parentless_pipeline_child_has_upper_targets(result, cluster_id);
            parentless_pipeline_export_or_extract_stopped_cluster(
                result,
                cluster_id,
                config,
                imported_source_cluster_ids_by_target[cluster_id],
                advance.touched_active_detector_sources,
                direct_result_superseded,
                has_upper_targets);
        }
        if (!runnable_cluster_ids.empty()) {
            prior_level_critical_events = level_max_critical_events;
            prior_level_stop_time = std::max(prior_level_stop_time, level_max_stop_time);
            have_prior_completed_level = true;
        }
        parentless_pipeline_queue_finished_level_imports(
            result, level, pending_import_sources_by_target);
    }
    result.profiling_stats.scheduler_wall_ns += elapsed_ns(scheduler_start);
    aggregate_parentless_pipeline_roots(graph, result, config, direct_result_superseded);
    finalize_cluster_step_parallel_accounting(result, config);
    return true;
}

}  // namespace

void prewarm_parentless_level_worker_pool(size_t worker_count) {
    if (worker_count > 1) {
        (void)parentless_level_worker_pool(worker_count);
    }
}

LockstepDecodeResult lockstep_hierarchical_decode(
    UserGraph& graph,
    const std::vector<uint64_t>& active_detectors,
    const LockstepSchedulerConfig& config) {
    const bool debug_decode_stages =
        std::getenv("PYMATCHING_DEBUG_TRACE_DECODE_STAGES") != nullptr &&
        std::string(std::getenv("PYMATCHING_DEBUG_TRACE_DECODE_STAGES")) == "1";
    if (debug_decode_stages) {
        std::cerr << "PYMATCHING DEBUG DECODE ENTER active=" << active_detectors.size() << std::endl;
    }
    if (config.edge_correlations) {
        throw std::invalid_argument("Lockstep hierarchical decoding does not yet support edge correlations.");
    }
    checked_num_workers(config);

    LockstepDecodeResult result;
    result.root_aggregate_result = ExtendedMatchingResult(graph.get_num_observables());

    if (config.enable_root_global_decode_fast_path) {
        const char* keep_hierarchy_env = std::getenv("PYMATCHING_NOSUBGRAPH_KEEP_HIERARCHY");
        const bool keep_hierarchy = keep_hierarchy_env != nullptr && std::string(keep_hierarchy_env) == "1";
        if (!keep_hierarchy) {
            // Fully lazy no-subgraph execution. If every worker receives the
            // original detector graph and the assigned active detector set, the
            // per-shot hierarchy only needs to name the active set.  In the
            // canonical default case there is one max/root worker, so avoid the
            // old influence/parent/subgraph materialization entirely.  The
            // benchmark still validates every observable prediction against the
            // ordinary global decoder.
            if (config.edge_correlations || graph.get_num_observables() > sizeof(pm::obs_int) * 8) {
                throw std::invalid_argument(
                    "lazy no-subgraph fast path currently supports only non-correlated <=64-observable decoding.");
            }
            ProcessingCluster root;
            root.id = 0;
            root.level = 0;
            root.active_detectors = active_detectors;
            root.influence_is_full_graph = true;
            root.parent_id = NO_PROCESSING_CLUSTER_PARENT;
            root.diameter_bound = 0;
            root.buffer_bound = 0;
            root.diameter = 0;
            result.hierarchy.clusters.push_back(std::move(root));
            result.hierarchy.cluster_ids_by_level.push_back({0});
            result.hierarchy.root_cluster_ids.push_back(0);
            result.profiling_stats.cluster_count = 1;
            result.profiling_stats.root_cluster_count = 1;
            result.profiling_stats.max_level = 0;
            result.profiling_stats.max_clusters_in_level = 1;
            result.profiling_stats.sum_cluster_active_detectors = active_detectors.size();
            result.profiling_stats.max_cluster_active_detectors = active_detectors.size();

            auto& mwpm = graph.get_mwpm();
            auto decode_start = steady_clock::now();
            auto compact = decode_detection_events_for_up_to_64_observables(
                mwpm, active_detectors, /*edge_correlations=*/false);
            auto decode_elapsed = elapsed_ns(decode_start);
            result.profiling_stats.single_root_decode_wall_ns = decode_elapsed;
            result.profiling_stats.cluster_step_wall_ns = decode_elapsed;
            result.profiling_stats.max_cluster_step_wall_ns = decode_elapsed;
            result.profiling_stats.max_level_cluster_step_wall_ns = decode_elapsed;
            result.profiling_stats.max_level_cluster_lifecycle_wall_ns = decode_elapsed;
            result.profiling_stats.max_level_cluster_id = 0;
            result.profiling_stats.ideal_worker_cluster_makespan_wall_ns = decode_elapsed;
            result.profiling_stats.ideal_cluster_worker_count_used =
                config.ideal_cluster_worker_count ? config.ideal_cluster_worker_count : config.num_workers;
            fill_bit_vector_from_obs_mask(
                compact.obs_mask,
                result.root_aggregate_result.obs_crossed.data(),
                graph.get_num_observables());
            result.root_aggregate_result.weight = compact.weight;
            return result;
        }
    }

    auto* graph_cache = graph.get_processing_cluster_graph_cache();
    const auto& matching_graph = graph.get_matching_graph_for_parallel_clustering();
    if (debug_decode_stages) {
        std::cerr << "PYMATCHING DEBUG DECODE BEFORE_RADIUS_PRECOMPUTE active=" << active_detectors.size() << std::endl;
    }
    ensure_processing_cluster_radius_neighbors_precomputed(matching_graph, *graph_cache, config.cluster_config);
    if (debug_decode_stages) {
        std::cerr << "PYMATCHING DEBUG DECODE AFTER_RADIUS_PRECOMPUTE active=" << active_detectors.size() << std::endl;
        std::cerr << "PYMATCHING DEBUG DECODE BEFORE_HIERARCHY active=" << active_detectors.size() << std::endl;
    }
    auto clustering_start = steady_clock::now();
    ProcessingClusterProfilingStats clustering_profiling_stats;
    if (config.prebuilt_processing_cluster_hierarchy != nullptr) {
        result.hierarchy = *config.prebuilt_processing_cluster_hierarchy;
    } else if (config.enable_growing_only_processing_clusters) {
        result.hierarchy = build_growing_only_processing_cluster_hierarchy(
            matching_graph, active_detectors, config.cluster_config, graph_cache, &clustering_profiling_stats);
    } else {
        result.hierarchy =
            build_processing_cluster_hierarchy(
                matching_graph, active_detectors, config.cluster_config, graph_cache, &clustering_profiling_stats);
        expand_parent_influence_to_cover_direct_children(result.hierarchy, graph.get_num_nodes());
    }
    result.profiling_stats.clustering_wall_ns += elapsed_ns(clustering_start);
    if (debug_decode_stages) {
        std::cerr << "PYMATCHING DEBUG DECODE AFTER_HIERARCHY active=" << active_detectors.size()
                  << " clusters=" << result.hierarchy.clusters.size()
                  << " roots=" << result.hierarchy.root_cluster_ids.size()
                  << " levels=" << result.hierarchy.cluster_ids_by_level.size()
                  << " clustering_ms=" << (static_cast<double>(result.profiling_stats.clustering_wall_ns) / 1e6)
                  << std::endl;
    }
    result.profiling_stats.clustering_component_construction_wall_ns +=
        clustering_profiling_stats.component_construction_wall_ns;
    result.profiling_stats.clustering_diameter_check_wall_ns += clustering_profiling_stats.diameter_check_wall_ns;
    result.profiling_stats.clustering_precomputable_distance_lookup_wall_ns +=
        clustering_profiling_stats.precomputable_distance_lookup_wall_ns;
    result.profiling_stats.clustering_influence_region_wall_ns +=
        clustering_profiling_stats.influence_region_wall_ns;
    result.profiling_stats.clustering_parent_assignment_wall_ns +=
        clustering_profiling_stats.parent_assignment_wall_ns;
    result.profiling_stats.growing_only_used_active_pair_table +=
        clustering_profiling_stats.growing_only_used_active_pair_table;
    result.profiling_stats.growing_only_used_sparse_frontier +=
        clustering_profiling_stats.growing_only_used_sparse_frontier;
    result.profiling_stats.growing_only_used_component_lookup +=
        clustering_profiling_stats.growing_only_used_component_lookup;
    result.profiling_stats.growing_only_forced_max_level_shortcut_count +=
        clustering_profiling_stats.growing_only_forced_max_level_shortcut_count;
    result.profiling_stats.growing_only_processed_frontier_events +=
        clustering_profiling_stats.growing_only_processed_frontier_events;
    result.profiling_stats.growing_only_generated_collision_events +=
        clustering_profiling_stats.growing_only_generated_collision_events;
    result.profiling_stats.growing_only_processed_collision_events +=
        clustering_profiling_stats.growing_only_processed_collision_events;
    result.profiling_stats.growing_only_exact_internal_pair_checks +=
        clustering_profiling_stats.growing_only_exact_internal_pair_checks;
    result.profiling_stats.growing_only_exact_external_pair_checks +=
        clustering_profiling_stats.growing_only_exact_external_pair_checks;
    result.profiling_stats.growing_only_exact_diameter_checks +=
        clustering_profiling_stats.growing_only_exact_diameter_checks;
    result.profiling_stats.growing_only_max_component_active_size = std::max<uint64_t>(
        result.profiling_stats.growing_only_max_component_active_size,
        clustering_profiling_stats.growing_only_max_component_active_size);
    result.profiling_stats.growing_only_max_residual_active_size = std::max<uint64_t>(
        result.profiling_stats.growing_only_max_residual_active_size,
        clustering_profiling_stats.growing_only_max_residual_active_size);
    result.profiling_stats.cluster_count = result.hierarchy.clusters.size();
    result.profiling_stats.root_cluster_count = result.hierarchy.root_cluster_ids.size();
    result.profiling_stats.max_level = find_maximum_level(result.hierarchy);
    result.profiling_stats.max_clusters_in_level = find_max_clusters_in_level(result.hierarchy);
    result.profiling_stats.sum_cluster_active_detectors = sum_cluster_active_detectors(result.hierarchy);
    result.profiling_stats.max_cluster_active_detectors = max_cluster_active_detectors(result.hierarchy);

    result.direct_child_cluster_ids.assign(result.hierarchy.clusters.size(), {});
    for (size_t cluster_id = 0; cluster_id < result.hierarchy.clusters.size(); cluster_id++) {
        auto parent_id = result.hierarchy.clusters[cluster_id].parent_id;
        if (parent_id != NO_PROCESSING_CLUSTER_PARENT) {
            result.direct_child_cluster_ids[parent_id].push_back(cluster_id);
        }
    }
    result.cluster_execution_states.assign(
        result.hierarchy.clusters.size(), LockstepClusterExecutionState::RUNNING);
    result.cluster_has_stopped_state.assign(result.hierarchy.clusters.size(), false);
    result.cluster_result_is_provisional.assign(result.hierarchy.clusters.size(), false);
    result.cluster_has_overlap_checkpoint.assign(result.hierarchy.clusters.size(), false);
    result.cluster_overlap_checkpoint_times.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_overlap_checkpoints.resize(result.hierarchy.clusters.size());
    result.cluster_overlap_checkpoint_extraction_detection_events.resize(result.hierarchy.clusters.size());
    result.cluster_used_overlap_checkpoint_recovery.assign(result.hierarchy.clusters.size(), false);
    result.cluster_overlap_snapshot_histories.resize(result.hierarchy.clusters.size());
    result.cluster_last_progress_time.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_extraction_detection_events.resize(result.hierarchy.clusters.size());
    result.cluster_extraction_detection_event_seen.resize(result.hierarchy.clusters.size());
    result.cluster_stopped_states.resize(result.hierarchy.clusters.size());
    result.cluster_direct_compact_results.assign(result.hierarchy.clusters.size(), MatchingResult());
    result.cluster_has_direct_compact_result.assign(result.hierarchy.clusters.size(), false);
    result.cluster_debug_compact_results.assign(result.hierarchy.clusters.size(), MatchingResult());
    result.cluster_has_debug_compact_result.assign(result.hierarchy.clusters.size(), false);
    result.cluster_direct_results.assign(
        result.hierarchy.clusters.size(), ExtendedMatchingResult(graph.get_num_observables()));
    result.cluster_has_direct_result.assign(result.hierarchy.clusters.size(), false);
    result.cluster_step_wall_ns.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_step_event_counts.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_critical_path_event_counts.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_has_critical_path_event_count.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_has_cutoff_algorithmic_time.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_cutoff_algorithmic_times.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_prior_level_critical_event_counts.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_overlap_event_bases.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_events_at_or_before_cutoff.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_events_after_cutoff.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_stop_algorithmic_times.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_critical_path_algorithmic_times.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_initial_setup_wall_ns.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_child_import_wall_ns.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_checkpoint_capture_wall_ns.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_checkpoint_restore_wall_ns.assign(result.hierarchy.clusters.size(), 0);
    result.cluster_mark_stop_wall_ns.assign(result.hierarchy.clusters.size(), 0);

    if (config.enable_root_global_decode_fast_path) {
        // No-subgraph execution model.  Under the invariant that a worker's
        // blossom never needs to grow outside its assigned root/max-level active
        // component, the worker can run directly on the immutable detector graph
        // and only needs the active detector list.  Avoid ClusterSubgraph
        // materialization, local/global id maps, boundary remapping, stopped-state
        // import/export, and scheduler bookkeeping.  The benchmark validates the
        // observable prediction against the ordinary global decoder.
        if (config.edge_correlations || graph.get_num_observables() > sizeof(pm::obs_int) * 8) {
            throw std::invalid_argument(
                "root-global no-subgraph fast path currently supports only non-correlated <=64-observable decoding.");
        }
        pm::MatchingResult compact_root_aggregate;
        uint64_t max_root_decode_ns = 0;
        for (auto root_cluster_id : result.hierarchy.root_cluster_ids) {
            if (root_cluster_id >= result.hierarchy.clusters.size()) {
                continue;
            }
            const auto& root_cluster = result.hierarchy.clusters[root_cluster_id];
            std::vector<uint64_t> root_detection_events;
            if (result.hierarchy.root_cluster_ids.size() == 1) {
                // The single root receives the original shot syndrome.  In the
                // old subgraph scheduler, child stop states would be imported
                // into this root before extraction; when subgraphs/imports are
                // intentionally removed, the equivalent no-subgraph input is the
                // full active-detector set.
                root_detection_events = active_detectors;
            } else {
                for (size_t cid = 0; cid < result.hierarchy.clusters.size(); cid++) {
                    size_t cursor = cid;
                    while (cursor != NO_PROCESSING_CLUSTER_PARENT && cursor < result.hierarchy.clusters.size() &&
                           cursor != root_cluster_id) {
                        cursor = result.hierarchy.clusters[cursor].parent_id;
                    }
                    if (cursor == root_cluster_id) {
                        const auto& a = result.hierarchy.clusters[cid].active_detectors;
                        root_detection_events.insert(root_detection_events.end(), a.begin(), a.end());
                    }
                }
                std::sort(root_detection_events.begin(), root_detection_events.end());
                root_detection_events.erase(
                    std::unique(root_detection_events.begin(), root_detection_events.end()),
                    root_detection_events.end());
            }
            auto& mwpm = graph.get_mwpm();
            auto root_decode_start = steady_clock::now();
            auto compact = decode_detection_events_for_up_to_64_observables(
                mwpm, root_detection_events, /*edge_correlations=*/false);
            auto root_decode_elapsed = elapsed_ns(root_decode_start);
            max_root_decode_ns = std::max(max_root_decode_ns, root_decode_elapsed);
            compact_root_aggregate += compact;
            if (root_cluster_id < result.cluster_step_wall_ns.size()) {
                result.cluster_step_wall_ns[root_cluster_id] += root_decode_elapsed;
            }
            if (root_cluster_id < result.cluster_extraction_detection_events.size()) {
                result.cluster_extraction_detection_events[root_cluster_id] = root_detection_events;
            }
            if (root_cluster.level == result.profiling_stats.max_level) {
                result.profiling_stats.max_level_cluster_step_wall_ns =
                    std::max(result.profiling_stats.max_level_cluster_step_wall_ns, root_decode_elapsed);
                if (root_decode_elapsed >= result.profiling_stats.max_level_cluster_lifecycle_wall_ns) {
                    result.profiling_stats.max_level_cluster_lifecycle_wall_ns = root_decode_elapsed;
                    result.profiling_stats.max_level_cluster_id = root_cluster_id;
                }
            }
            result.cluster_execution_states[root_cluster_id] = LockstepClusterExecutionState::QUIESCENT_STOPPED;
            result.cluster_has_stopped_state[root_cluster_id] = true;
            if (root_cluster_id < result.cluster_has_direct_compact_result.size()) {
                result.cluster_direct_compact_results[root_cluster_id] = compact;
                result.cluster_has_direct_compact_result[root_cluster_id] = true;
            }
        }
        result.profiling_stats.single_root_decode_wall_ns += max_root_decode_ns;
        result.profiling_stats.cluster_step_wall_ns += max_root_decode_ns;
        fill_bit_vector_from_obs_mask(
            compact_root_aggregate.obs_mask,
            result.root_aggregate_result.obs_crossed.data(),
            graph.get_num_observables());
        result.root_aggregate_result.weight = compact_root_aggregate.weight;
        finalize_cluster_step_parallel_accounting(result, config);
        return result;
    }

    if (debug_decode_stages) {
        std::cerr << "PYMATCHING DEBUG DECODE BEFORE_SUBGRAPH_BUILD active=" << active_detectors.size()
                  << " clusters=" << result.hierarchy.clusters.size() << std::endl;
    }
    auto subgraph_build_start = steady_clock::now();
    std::vector<uint64_t> cluster_build_wall_ns;
    result.cluster_subgraphs =
        build_cluster_subgraphs_parallel(graph, result.hierarchy.clusters, config, &cluster_build_wall_ns);
    auto subgraph_build_elapsed = elapsed_ns(subgraph_build_start);
    if (debug_decode_stages) {
        std::cerr << "PYMATCHING DEBUG DECODE AFTER_SUBGRAPH_BUILD active=" << active_detectors.size()
                  << " subgraph_build_ms=" << (static_cast<double>(subgraph_build_elapsed) / 1e6)
                  << " clusters=" << result.hierarchy.clusters.size() << std::endl;
    }
    if (config.enable_full_graph_worker_subgraphs) {
        uint64_t build_total = 0;
        for (size_t k = 0; k < cluster_build_wall_ns.size(); k++) {
            build_total += cluster_build_wall_ns[k];
        }
        result.profiling_stats.worker_graph_distribution_wall_ns += build_total;
        if (!config.exclude_full_graph_worker_distribution_from_runtime) {
            // Diagnostic accounting mode only: charge worker-local full graph
            // materialization to the cluster lifecycle.  The default policy
            // treats these detector graph copies as pre-distributed to workers
            // before online decoding starts, so this block is normally skipped.
            for (size_t k = 0; k < cluster_build_wall_ns.size() && k < result.cluster_initial_setup_wall_ns.size(); k++) {
                result.cluster_initial_setup_wall_ns[k] += cluster_build_wall_ns[k];
            }
            result.profiling_stats.initial_mwpm_setup_wall_ns += build_total;
        }
    } else {
        result.profiling_stats.subgraph_build_wall_ns += subgraph_build_elapsed;
    }

    if (config.enable_parentless_level_pipeline) {
        if (debug_decode_stages) {
            std::cerr << "PYMATCHING DEBUG DECODE BEFORE_PARENTLESS_PIPELINE active=" << active_detectors.size()
                      << " clusters=" << result.hierarchy.clusters.size() << std::endl;
        }
        try {
            if (run_parentless_level_pipeline(graph, result, config)) {
                if (debug_decode_stages) {
                    std::cerr << "PYMATCHING DEBUG DECODE AFTER_PARENTLESS_PIPELINE active=" << active_detectors.size()
                              << " clusters=" << result.hierarchy.clusters.size() << std::endl;
                }
                return result;
            }
        } catch (const std::invalid_argument&) {
            if (!config.allow_global_fallback) {
                throw;
            }
            result.used_global_fallback = true;
            result.root_aggregate_result = decode_full_graph(graph, active_detectors, config.edge_correlations);
            return result;
        }
    }

    if (result.hierarchy.clusters.size() == 1 && result.hierarchy.root_cluster_ids.size() == 1) {
        auto& root_subgraph = result.cluster_subgraphs[0];
        auto& root_mwpm = cluster_mwpm(root_subgraph, config.edge_correlations);
        result.cluster_extraction_detection_events[0] = root_subgraph.local_active_detectors;
        try {
            auto root_decode_start = steady_clock::now();
            process_timeline_until_completion(root_mwpm, root_subgraph.local_active_detectors);
            result.profiling_stats.single_root_decode_wall_ns += elapsed_ns(root_decode_start);
        } catch (const std::invalid_argument&) {
            if (!config.allow_global_fallback) {
                throw;
            }
            result.used_global_fallback = true;
            result.root_aggregate_result = decode_full_graph(graph, active_detectors, config.edge_correlations);
            return result;
        }
        result.cluster_execution_states[0] = LockstepClusterExecutionState::QUIESCENT_STOPPED;
        result.cluster_has_stopped_state[0] = true;
        auto root_extract_start = steady_clock::now();
        if (root_mwpm.flooder.graph.num_observables <= sizeof(pm::obs_int) * 8) {
            result.cluster_direct_compact_results[0] =
                extract_compact_result_from_current_mwpm_state(root_mwpm, result.cluster_extraction_detection_events[0]);
            result.cluster_has_direct_compact_result[0] = true;
            fill_bit_vector_from_obs_mask(
                result.cluster_direct_compact_results[0].obs_mask,
                result.root_aggregate_result.obs_crossed.data(),
                graph.get_num_observables());
            result.root_aggregate_result.weight = result.cluster_direct_compact_results[0].weight;
            if (root_subgraph.use_shared_full_graph_mwpm) {
                root_mwpm.reset();
            }
        } else {
            result.cluster_direct_results[0] =
                extract_result_from_current_mwpm_state(root_mwpm, result.cluster_extraction_detection_events[0]);
            result.cluster_has_direct_result[0] = true;
            result.root_aggregate_result = result.cluster_direct_results[0];
            if (root_subgraph.use_shared_full_graph_mwpm) {
                root_mwpm.reset();
            }
        }
        result.profiling_stats.root_extraction_wall_ns += elapsed_ns(root_extract_start);
        return result;
    }

    auto initial_setup_start = steady_clock::now();
    for (size_t cluster_id = 0; cluster_id < result.cluster_subgraphs.size(); cluster_id++) {
        auto cluster_init_start = steady_clock::now();
        auto& subgraph = result.cluster_subgraphs[cluster_id];
        result.cluster_extraction_detection_events[cluster_id] = subgraph.local_active_detectors;
        rebuild_extraction_seen_for_cluster(result, cluster_id);
        auto& mwpm = cluster_mwpm(subgraph, config.edge_correlations);
        if (subgraph.persistent_mwpm != nullptr && !subgraph.persistent_mwpm_acquired_clean) {
            mwpm.reset();
        }
        initialize_mwpm_for_detection_events(mwpm, subgraph.local_active_detectors);
        subgraph.persistent_mwpm_acquired_clean = false;
        if (cluster_id < result.cluster_initial_setup_wall_ns.size()) {
            result.cluster_initial_setup_wall_ns[cluster_id] += elapsed_ns(cluster_init_start);
        }
        maybe_mark_cluster_stopped(result, cluster_id, config);
    }
    result.profiling_stats.initial_mwpm_setup_wall_ns += elapsed_ns(initial_setup_start);
    for (size_t cluster_id = 0; cluster_id < result.cluster_subgraphs.size(); cluster_id++) {
        if (result.cluster_execution_states[cluster_id] == LockstepClusterExecutionState::RUNNING ||
            result.cluster_execution_states[cluster_id] == LockstepClusterExecutionState::QUIESCENT_STOPPED ||
            result.cluster_execution_states[cluster_id] == LockstepClusterExecutionState::DRAINED_WITH_LIVE_ALT_TREE) {
            maybe_capture_overlap_checkpoint_for_parent(result, config, cluster_id);
        }
    }
    try {
        mark_stopped_running_clusters_and_import_until_stable(graph, result, config);
    } catch (const std::invalid_argument&) {
        if (!config.allow_global_fallback) {
            throw;
        }
        result.used_global_fallback = true;
        result.root_aggregate_result = decode_full_graph(graph, active_detectors, config.edge_correlations);
        return result;
    }

    auto scheduler_start = steady_clock::now();
    while (true) {
        try {
            mark_stopped_running_clusters_and_import_until_stable(graph, result, config);
        } catch (const std::invalid_argument&) {
            if (!config.allow_global_fallback) {
                throw;
            }
            result.used_global_fallback = true;
            result.root_aggregate_result = decode_full_graph(graph, active_detectors, config.edge_correlations);
            return result;
        }
        auto ready_scan = scan_clusters_at_next_global_time(result, config);
        if (!ready_scan.found) {
            break;
        }
        auto next_global_time = ready_scan.time;
        while (true) {
            for (size_t cluster_id = 0; cluster_id < result.cluster_subgraphs.size(); cluster_id++) {
                if (result.cluster_execution_states[cluster_id] != LockstepClusterExecutionState::RUNNING) {
                    continue;
                }
                auto& subgraph = result.cluster_subgraphs[cluster_id];
                auto& mwpm = cluster_mwpm(subgraph, config.edge_correlations);
                if (!mwpm.flooder.has_valid_tentative_events()) {
                    maybe_mark_cluster_stopped(result, cluster_id, config);
                    continue;
                }
            }
            auto cluster_ids_to_step = collect_clusters_to_step_at_time(result, config, next_global_time);
            if (!cluster_ids_to_step.empty()) {
                result.profiling_stats.scheduler_batch_count++;
                result.profiling_stats.scheduler_cluster_step_count += cluster_ids_to_step.size();
            }
            for (auto cluster_id : cluster_ids_to_step) {
                maybe_capture_overlap_checkpoint_for_parent(result, config, cluster_id);
            }
            auto step_outcomes = step_clusters_at_time(result, config, cluster_ids_to_step, next_global_time);
            for (const auto& outcome : step_outcomes) {
                result.cluster_last_progress_time[outcome.cluster_id] = outcome.progress_time;
                maybe_mark_cluster_stopped(result, outcome.cluster_id, config);
                capture_overlap_history_snapshot_for_cluster(result, config, outcome.cluster_id);
            }
            try {
                mark_stopped_running_clusters_and_import_until_stable(graph, result, config);
            } catch (const std::invalid_argument&) {
                if (!config.allow_global_fallback) {
                    throw;
                }
                result.used_global_fallback = true;
                result.root_aggregate_result = decode_full_graph(graph, active_detectors, config.edge_correlations);
                return result;
            }
            if (!has_running_cluster_with_scheduled_work_at_time(result, config, next_global_time)) {
                break;
            }
        }
    }
    try {
        mark_stopped_running_clusters_and_import_until_stable(graph, result, config);
    } catch (const std::invalid_argument&) {
        if (!config.allow_global_fallback) {
            throw;
        }
        result.used_global_fallback = true;
        result.root_aggregate_result = decode_full_graph(graph, active_detectors, config.edge_correlations);
        return result;
    }
    result.profiling_stats.scheduler_wall_ns += elapsed_ns(scheduler_start);
    finalize_cluster_step_parallel_accounting(result, config);

    bool has_non_root_provisional_result = false;
    bool hierarchy_requires_deep_import = false;
    for (size_t cluster_id = 0; cluster_id < result.cluster_subgraphs.size(); cluster_id++) {
        if (result.cluster_execution_states[cluster_id] == LockstepClusterExecutionState::RUNNING) {
            maybe_mark_cluster_stopped(result, cluster_id, config);
        }
        if (result.cluster_result_is_provisional[cluster_id]) {
            has_non_root_provisional_result = true;
        }
        if (result.hierarchy.clusters[cluster_id].parent_id != NO_PROCESSING_CLUSTER_PARENT &&
            result.cluster_execution_states[cluster_id] != LockstepClusterExecutionState::IMPORTED_INTO_PARENT) {
            hierarchy_requires_deep_import = true;
        }
    }

    if (result.encountered_non_quiescent_cluster || has_non_root_provisional_result || hierarchy_requires_deep_import) {
        if (!config.allow_global_fallback) {
            size_t running=0, quiescent=0, drained=0, imported=0, provisional=0, nonroot_not_imported=0;
            for (size_t cid = 0; cid < result.cluster_execution_states.size(); cid++) {
                auto st = result.cluster_execution_states[cid];
                if (st == LockstepClusterExecutionState::RUNNING) running++;
                else if (st == LockstepClusterExecutionState::QUIESCENT_STOPPED) quiescent++;
                else if (st == LockstepClusterExecutionState::DRAINED_WITH_LIVE_ALT_TREE) drained++;
                else if (st == LockstepClusterExecutionState::IMPORTED_INTO_PARENT) imported++;
                if (result.cluster_result_is_provisional[cid]) provisional++;
                if (result.hierarchy.clusters[cid].parent_id != NO_PROCESSING_CLUSTER_PARENT && st != LockstepClusterExecutionState::IMPORTED_INTO_PARENT) nonroot_not_imported++;
            }
            throw std::invalid_argument(
                "Lockstep hierarchical decoder needs state import for a nested or non-quiescent cluster, but global fallback is disabled. counts running=" + std::to_string(running) +
                " quiescent=" + std::to_string(quiescent) + " drained=" + std::to_string(drained) +
                " imported=" + std::to_string(imported) + " provisional=" + std::to_string(provisional) +
                " nonroot_not_imported=" + std::to_string(nonroot_not_imported));
        }
        result.used_global_fallback = true;
        result.root_aggregate_result = decode_full_graph(graph, active_detectors, config.edge_correlations);
        return result;
    }

    const bool can_use_compact_root_aggregation = graph.get_num_observables() <= sizeof(pm::obs_int) * 8;
    pm::MatchingResult compact_root_aggregate;

    for (auto root_cluster_id : result.hierarchy.root_cluster_ids) {
        if (!result.cluster_has_stopped_state[root_cluster_id]) {
            if (!config.allow_global_fallback) {
                throw std::invalid_argument(
                    "Lockstep hierarchical decoder could not stop every root cluster, but global fallback is disabled. root=" + std::to_string(root_cluster_id) +
                    " state=" + std::to_string(static_cast<int>(result.cluster_execution_states[root_cluster_id])) +
                    " has_state=" + std::to_string(result.cluster_has_stopped_state[root_cluster_id]) +
                    " scheduled=" + std::to_string(cluster_has_scheduled_work(result.cluster_subgraphs[root_cluster_id], config.edge_correlations)));
            }
            result.used_global_fallback = true;
            result.root_aggregate_result = decode_full_graph(graph, active_detectors, config.edge_correlations);
            return result;
        }

        if (can_use_compact_root_aggregation) {
            if (root_cluster_id < result.cluster_has_direct_compact_result.size() &&
                result.cluster_has_direct_compact_result[root_cluster_id]) {
                compact_root_aggregate += result.cluster_direct_compact_results[root_cluster_id];
                continue;
            }

            auto& root_subgraph = result.cluster_subgraphs[root_cluster_id];
            auto& root_mwpm = cluster_mwpm(root_subgraph, config.edge_correlations);
            auto root_extract_start = steady_clock::now();
            compact_root_aggregate += extract_compact_result_from_current_mwpm_state(
                root_mwpm, result.cluster_extraction_detection_events[root_cluster_id]);
            result.profiling_stats.root_extraction_wall_ns += elapsed_ns(root_extract_start);
            if (root_subgraph.persistent_mwpm != nullptr) {
                release_persistent_full_graph_worker_mwpm(root_subgraph.persistent_mwpm);
                root_subgraph.persistent_mwpm = nullptr;
            } else if (root_subgraph.use_shared_full_graph_mwpm) {
                root_mwpm.reset();
            }
            continue;
        }

        if (root_cluster_id < result.cluster_has_direct_result.size() &&
            result.cluster_has_direct_result[root_cluster_id]) {
            result.root_aggregate_result += result.cluster_direct_results[root_cluster_id];
            continue;
        }

        auto& root_subgraph = result.cluster_subgraphs[root_cluster_id];
        auto& root_mwpm = cluster_mwpm(root_subgraph, config.edge_correlations);
        auto root_extract_start = steady_clock::now();
        result.root_aggregate_result += extract_result_from_current_mwpm_state(
            root_mwpm, result.cluster_extraction_detection_events[root_cluster_id]);
        result.profiling_stats.root_extraction_wall_ns += elapsed_ns(root_extract_start);
        if (root_subgraph.persistent_mwpm != nullptr) {
            release_persistent_full_graph_worker_mwpm(root_subgraph.persistent_mwpm);
            root_subgraph.persistent_mwpm = nullptr;
        } else if (root_subgraph.use_shared_full_graph_mwpm) {
            root_mwpm.reset();
        }
    }

    if (can_use_compact_root_aggregation) {
        fill_bit_vector_from_obs_mask(
            compact_root_aggregate.obs_mask,
            result.root_aggregate_result.obs_crossed.data(),
            graph.get_num_observables());
        result.root_aggregate_result.weight = compact_root_aggregate.weight;
    }

    return result;
}

}  // namespace pm
