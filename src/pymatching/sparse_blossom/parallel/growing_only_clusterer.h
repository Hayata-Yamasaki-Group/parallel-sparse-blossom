// Copyright 2026
// Experimental sparse-blossom-like growing-only clustering path.
//
// This file intentionally does not participate in MWPM matching.  It reuses the
// detector graph metric and the same ProcessingClusterConfig level bounds, but
// replaces blossom/matching state updates with a small UF over active detector
// regions.  The purpose is online preprocessing experiments only.

#ifndef PYMATCHING_PARALLEL_GROWING_ONLY_CLUSTERER_H
#define PYMATCHING_PARALLEL_GROWING_ONLY_CLUSTERER_H

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "pymatching/sparse_blossom/flooder/graph.h"
#include "pymatching/sparse_blossom/ints.h"
#include "pymatching/sparse_blossom/parallel/processing_cluster.h"

namespace pm {

constexpr cumulative_time_int GROWING_ONLY_INF_DISTANCE = std::numeric_limits<cumulative_time_int>::max() / 4;

struct GrowingOnlyClusterMeta {
    uint32_t uf_parent = 0;
    uint32_t size = 1;

    // Max time of the merge events needed to form this connected component.
    cumulative_time_int ready_time = 0;

    // Max same-UF active-pair collision time observed so far.  This is updated
    // for same-component collisions after the component was already connected;
    // it is not limited to the merge edge set.
    cumulative_time_int max_internal_pair_collision_time = 0;

    // Exact detector-graph distance corresponding to the maximum internal pair
    // collision. Existing ProcessingCluster diameter checks use this distance
    // unit directly. For normal even-weight sparse-blossom graphs this equals
    // 2 * max_internal_pair_collision_time.
    cumulative_time_int max_internal_pair_distance = 0;

    // Nearest way for an odd-parity component to stop at the graph boundary.
    // Existing ProcessingCluster accepts an odd active-detector cluster only
    // when this distance is <= diameter_bound.
    cumulative_time_int nearest_boundary_match_distance = GROWING_ONLY_INF_DISTANCE;

    // First collision with an active detector outside the connected component,
    // observed after ready_time.  INF_DISTANCE means not yet observed.
    cumulative_time_int first_external_collision_time = GROWING_ONLY_INF_DISTANCE;

    bool accepted = false;
    bool removed_or_hidden = false;
};

struct GrowingOnlyAcceptedCluster {
    std::vector<uint64_t> active_detectors;
    cumulative_time_int diameter = 0;
    cumulative_time_int buffer = 0;
    cumulative_time_int ready_time = 0;
    cumulative_time_int max_internal_pair_collision_time = 0;
    cumulative_time_int max_internal_pair_distance = 0;
    cumulative_time_int first_external_collision_time = GROWING_ONLY_INF_DISTANCE;
    cumulative_time_int nearest_boundary_match_distance = GROWING_ONLY_INF_DISTANCE;
    bool parity_can_stop_locally = false;
    bool forced_by_max_level = false;
    size_t internal_collisions_recorded = 0;
    size_t external_collisions_recorded = 0;
};

struct GrowingOnlyClustererStats {
    uint64_t event_generation_wall_ns = 0;
    uint64_t event_processing_wall_ns = 0;
    size_t generated_pair_collision_events = 0;
    size_t processed_pair_collision_events = 0;
    size_t internal_collisions_recorded = 0;
    size_t external_collisions_recorded = 0;
    size_t union_count = 0;
    size_t removed_or_hidden_clusters = 0;
    size_t parity_rejected_clusters = 0;
    size_t diameter_rejected_clusters = 0;
    size_t max_active_rejected_clusters = 0;
    size_t adaptive_deferred_clusters = 0;
    size_t adaptive_deferred_active_detectors = 0;
    size_t forced_by_max_level_clusters = 0;
    size_t forced_max_level_shortcut_count = 0;
    size_t sparse_frontier_events = 0;
    size_t generated_collision_events = 0;
    size_t processed_frontier_events = 0;
    size_t processed_collision_events = 0;
    size_t exact_internal_pair_checks = 0;
    size_t exact_external_pair_checks = 0;
    size_t exact_diameter_checks = 0;
    size_t max_component_active_size = 0;
    size_t max_residual_active_size = 0;
    bool used_sparse_frontier = false;
    bool used_component_lookup = false;
    bool used_active_pair_table = false;
};

struct GrowingOnlyClustererResult {
    size_t level = 1;
    ProcessingClusterBounds bounds{};
    std::vector<GrowingOnlyAcceptedCluster> accepted_clusters;
    std::vector<uint64_t> residual_active_detectors;
    GrowingOnlyClustererStats stats;

    size_t residual_active_detector_count() const {
        return residual_active_detectors.size();
    }
};

class GrowingOnlyClusterer {
   public:
    GrowingOnlyClusterer(
        const MatchingGraph& graph,
        std::vector<uint64_t> active_detectors,
        size_t level,
        const ProcessingClusterConfig& config,
        const ProcessingClusterGraphCache* graph_cache = nullptr);

    GrowingOnlyClustererResult run();

   private:
    struct PairCollisionEvent {
        cumulative_time_int time = 0;
        cumulative_time_int distance = 0;
        uint32_t a = 0;
        uint32_t b = 0;
        bool operator<(const PairCollisionEvent& other) const;
    };

    const MatchingGraph& graph;
    std::vector<uint64_t> active_detectors;
    size_t level;
    ProcessingClusterConfig config;
    const ProcessingClusterGraphCache* graph_cache;
    ProcessingClusterBounds bounds;

    std::vector<GrowingOnlyClusterMeta> meta;
    std::vector<std::vector<uint32_t>> members_by_root;
    std::vector<size_t> internal_collision_count_by_root;
    std::vector<size_t> external_collision_count_by_root;
    GrowingOnlyClustererStats stats;

    std::vector<PairCollisionEvent> generate_pair_collision_events();
    std::vector<PairCollisionEvent> generate_sparse_frontier_collision_events();
    void run_sparse_frontier_streaming_unions();
    std::vector<cumulative_time_int> dijkstra_from_source(uint64_t source) const;
    cumulative_time_int distance_between_active_detectors(size_t a, size_t b) const;

    uint32_t find(uint32_t k);
    uint32_t unite(uint32_t a, uint32_t b, cumulative_time_int event_time, cumulative_time_int event_distance);
    void note_internal_collision(uint32_t root, cumulative_time_int event_time, cumulative_time_int event_distance);
    void note_external_collision(uint32_t root, cumulative_time_int event_time);
    cumulative_time_int nearest_boundary_match_distance_for_detector(uint64_t detector) const;
    cumulative_time_int farthest_boundary_match_distance_for_root(uint32_t root);
    bool parity_can_stop_locally(uint32_t root);
    bool force_at_max_level() const;
    size_t max_active_detectors_for_this_level() const;
    void evaluate_final_components(GrowingOnlyClustererResult& result);
    void emit_accepted_cluster(uint32_t root, GrowingOnlyClustererResult& result, bool forced_by_max_level);
    void collect_residual(GrowingOnlyClustererResult& result);
    bool try_run_no_union_singleton_fast_path(GrowingOnlyClustererResult& result);
    bool run_sparse_frontier_path(GrowingOnlyClustererResult& result);
    bool forced_phi_buffer_certificate_allows_component_split() const;
    bool run_forced_max_level_single_root_certificate_fallback(GrowingOnlyClustererResult& result);
    bool run_forced_max_level_direct_component_lookup_path(GrowingOnlyClustererResult& result);
    bool run_precomputed_component_lookup_path(GrowingOnlyClustererResult& result);
    bool run_active_pair_table_component_path(GrowingOnlyClustererResult& result);
    bool use_sparse_frontier_path() const;
    bool use_reference_active_pair_table_path() const;
    bool skip_forced_max_exact_metadata() const;
    bool compute_exact_external_metadata() const;
    void compute_exact_component_pair_metadata();
    void compute_first_external_collision_times();
};

GrowingOnlyClustererResult run_growing_only_clustering(
    const MatchingGraph& graph,
    std::vector<uint64_t> active_detectors,
    size_t level,
    const ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache = nullptr);

// Mutable-config wrapper used by adaptive scheduling.  If
// config.adaptive_cluster_bounds is true, this tries candidate (d_k,b_k) pairs
// for the current residual active detectors, writes the selected bounds into
// config.explicit_bounds_by_level[level-1], and returns the clustering result
// for those selected bounds.  If adaptive mode is disabled, it is equivalent to
// run_growing_only_clustering.
GrowingOnlyClustererResult run_growing_only_clustering_with_adaptive_bounds(
    const MatchingGraph& graph,
    std::vector<uint64_t> active_detectors,
    size_t level,
    ProcessingClusterConfig& config,
    const ProcessingClusterGraphCache* graph_cache = nullptr);

}  // namespace pm

#endif  // PYMATCHING_PARALLEL_GROWING_ONLY_CLUSTERER_H
