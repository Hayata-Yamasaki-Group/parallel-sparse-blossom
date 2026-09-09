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

#ifndef PYMATCHING_PARALLEL_PROCESSING_CLUSTER_H
#define PYMATCHING_PARALLEL_PROCESSING_CLUSTER_H

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <vector>

#include "pymatching/sparse_blossom/flooder/graph.h"
#include "pymatching/sparse_blossom/ints.h"

namespace pm {

constexpr size_t NO_PROCESSING_CLUSTER_PARENT = std::numeric_limits<size_t>::max();

struct ProcessingClusterBounds {
    cumulative_time_int diameter_bound = 0;
    cumulative_time_int buffer_bound = 0;

    // Optional exact schedule values used by experimental policies whose
    // mathematical d_k or b_k may exceed the safe runtime integer range.
    // Runtime code uses the clamped integer fields above; certificate checks
    // use these long-double schedule fields when present.
    bool has_long_double_schedule_bounds = false;
    long double schedule_diameter_bound = 0;
    long double schedule_buffer_bound = 0;
    bool runtime_bounds_clamped = false;

    // Optional certificate metadata for explicit phi-derived bounds.  The
    // actual buffer must satisfy all three lower bounds below before a
    // max-level forced split may emit multiple disjoint component roots.
    bool has_phi_buffer_certificate = false;
    double phi = 0;
    cumulative_time_int phi_required_buffer_bound = 0;
    // Optional lower bound used by the experimental phi-floor budget policy.
    // A value of 0 means the normal certificate did not use an extra budget.
    cumulative_time_int phi_budget_required_buffer_bound = 0;
    cumulative_time_int ratio_required_buffer_bound = 0;
    cumulative_time_int min_required_buffer_bound = 0;
};

struct AdaptiveEventTiming {
    uint64_t raw_event_count = 0;
    uint64_t events_after_cutoff = 0;
    cumulative_time_int stop_algorithmic_time = 0;
};

struct ProcessingClusterConfig {
    double beta = 72;
    double gamma = 291;
    double lambda = 320;
    size_t max_level = 64;
    bool force_cluster_at_max_level = true;

    // Experimental single-root completion policy.  Normal lower levels still
    // build as many locally-stoppable clusters as possible.  If parity-aware
    // construction leaves residual detection events at the effective maximum
    // level, append one forced full-graph final root at level max_level+1
    // instead of throwing.  This keeps the root MWPM state single-owner while
    // preserving the lower-level clustering benefits.
    bool force_single_final_root = true;

    // Optional explicit one-based clustering bounds.  Entry 0 is level 1,
    // entry 1 is level 2, and so on.  When this vector is non-empty, these
    // bounds are used instead of the paper-formula beta/gamma/lambda bounds.
    // The effective maximum level is min(max_level, explicit_bounds_by_level.size()).
    std::vector<ProcessingClusterBounds> explicit_bounds_by_level;

    // Optional one-based active-detector count caps.  Entry 0 applies to level 1,
    // entry 1 to level 2, and so on.  A value of 0 means no count cap for that
    // level.  When non-empty, a non-forced cluster candidate is accepted only
    // when the whole buffer-connected component contains at most this many
    // active detectors.  This is useful for experimental schedules such as
    // level k targeting 2k detectors.
    std::vector<size_t> max_active_detectors_by_level;

    // When true, validate explicit_bounds_by_level against the monotone
    // clustering preconditions and the stopping lemma inequality.  It is false
    // by default so old benchmarking parameter sweeps remain usable.
    bool enforce_stopping_lemma_bounds = false;
    // Used only by the root-direct benchmarking path. It lets the benchmark
    // jump directly to the empirically relevant maximum/root level for a fixed
    // parameter set instead of constructing all lower-level child clusters.
    size_t root_direct_min_level = 1;
    // Used only by the root-direct benchmarking path after global-output
    // verification. Skipping this check removes an O(component^2) pass.
    bool root_direct_skip_diameter_check = false;

    // Experimental adaptive bound mode.  When enabled, each shot/level selects
    // an integer (d_k,b_k) before clustering by trying candidate bounds on the
    // current residual active detectors and choosing the pair that maximizes the
    // number of accepted clusters under the existing clustering predicates.
    bool adaptive_cluster_bounds = false;
    size_t adaptive_candidate_limit = 64;
    bool adaptive_trace_bounds = false;
    // Strategy for adaptive bound selection.  The legacy "trial" strategy
    // samples a small d/b grid; the "breakpoint" strategy only tests b
    // values where the residual active-detector threshold components can
    // change, then derives the largest certificate-compatible d from b.
    bool adaptive_breakpoint_strategy = false;

    // d_k,b_k constraints for adaptive mode.  The environment reader enables
    // the legacy-gap and phi-certified schedule by default for correctness;
    // explicit DISABLE_* environment variables can still be used for ablation.
    bool adaptive_enforce_monotone_bounds = false;
    bool adaptive_enforce_legacy_gap_bound = false;
    bool adaptive_enforce_phi_buffer_certificate = false;
    double adaptive_min_buffer_to_diameter_ratio = 0.0;
    bool adaptive_min_buffer_ratio_first_level_only = false;
    double adaptive_phi_floor = 0.0;

    // Adaptive objective.  The historical objective maximizes the number of
    // accepted clusters at the current level.  The max-active objective instead
    // minimizes the largest accepted-cluster active-detector count, optionally
    // with a small lookahead over future adaptive levels.
    bool adaptive_objective_minimize_max_active = false;
    // Experimental objective that evaluates candidates by the actual MWPM event
    // count of their accepted active-detector sets when the benchmark supplies
    // adaptive_event_count_callback.  Without the callback it falls back to a
    // deterministic proxy based on active-detector count.
    bool adaptive_objective_minimize_max_events = false;
    // Objective aligned with the reported ideal parallel event metric. Candidate
    // rollouts minimize the all-level critical event proxy over the current
    // level and the lookahead levels: the maximum level-critical event count in
    // the rollout, with a final residual/root event count when the rollout
    // reaches the effective maximum level. This deliberately does not penalize
    // the mere existence of extra nonempty levels.
    bool adaptive_objective_minimize_ideal_events = false;
    // Historical weighting parameter. It is ignored by the exact ideal-events
    // objective because weighting would optimize a proxy different from the
    // reported ideal event-count metric.
    double adaptive_ideal_events_weight_power = 0.0;
    std::function<uint64_t(const std::vector<uint64_t>&)> adaptive_event_count_callback;
    // Exact scheduler-derived objective support.  When supplied, candidate
    // rollouts can evaluate the same level-by-level ideal event metric used by
    // the final benchmark output: run sparse blossom to completion, count only
    // notifications after the previous completed level's stop time, add the
    // level maximum, and update the stop-time cutoff.
    std::function<AdaptiveEventTiming(const std::vector<uint64_t>&, cumulative_time_int)> adaptive_event_timing_callback;
    uint64_t adaptive_exact_ideal_prior_critical_events = 0;
    cumulative_time_int adaptive_exact_ideal_prior_stop_time = 0;
    bool adaptive_exact_ideal_have_prior_completed_level = false;
    size_t adaptive_lookahead_levels = 0;

    // Experimental selective deferral.  When enabled in adaptive mode, an
    // accepted non-forced cluster whose actual event count or active-detector
    // count is above the configured threshold is moved back to the residual set
    // instead of being fixed at the current level.  This preserves the existing
    // correctness predicates by only delaying decisions to higher levels.
    bool adaptive_defer_heavy_clusters = false;
    size_t adaptive_defer_max_active = 0;
    uint64_t adaptive_defer_max_events = 0;

    // Experimental level-1 bound re-optimization.  Instead of accepting extra
    // short pairs after fixed level-1 clustering, try multiple candidate d_1
    // values, recompute the corresponding certified b_1, rerun level-1 from
    // scratch for each candidate, and choose by the adaptive objective plus
    // lookahead.  This makes patterns such as o-x-o and o-x-x-o part of the
    // actual level-1 certificate instead of a post-hoc heuristic.
    bool adaptive_level1_reoptimize_bounds = false;
    double adaptive_level1_max_diameter_ratio = 3.0;
    cumulative_time_int adaptive_level1_max_diameter = 0;
    size_t adaptive_level1_candidate_limit = 64;

    // Experimental level-1 short-pair completion.  After the ordinary level-1
    // clustering step, greedily accepts isolated residual pairs of active
    // detectors whose detector-graph distance is at most the configured bound.
    // This targets local patterns such as o-x-o and o-x-x-o without globally
    // increasing the level-1 buffer threshold and accidentally merging long
    // chains into large components.
    bool level1_short_pair_acceptance = false;
    cumulative_time_int level1_short_pair_max_distance = 0;
    double level1_short_pair_distance_ratio = 0.0;
    cumulative_time_int level1_short_pair_external_guard = 0;
    bool level1_short_pair_external_guard_explicit = false;
    bool level1_short_pair_require_mutual_nearest = true;
};

struct ProcessingCluster {
    size_t id;
    size_t level;
    std::vector<uint64_t> active_detectors;
    std::vector<size_t> influence_vertices;
    // True when influence_vertices is the full identity vertex set [0, num_nodes).
    // Keeping this as a marker avoids per-cluster O(|V|) vector materialization
    // for full-graph influence clusters; callers that need the explicit vector
    // can expand it from the graph size.
    bool influence_is_full_graph = false;
    std::vector<size_t> boundary_endpoint_vertices;
    size_t parent_id = NO_PROCESSING_CLUSTER_PARENT;
    cumulative_time_int diameter_bound;
    cumulative_time_int buffer_bound;
    cumulative_time_int diameter;
    cumulative_time_int nearest_boundary_match_distance = std::numeric_limits<cumulative_time_int>::max();
    cumulative_time_int nearest_external_detector_distance = std::numeric_limits<cumulative_time_int>::max();
    bool forced_by_max_level = false;
};

struct ProcessingClusterHierarchy {
    std::vector<ProcessingCluster> clusters;
    std::vector<std::vector<size_t>> cluster_ids_by_level;
    std::vector<size_t> root_cluster_ids;
};

struct ProcessingClusterProfilingStats {
    uint64_t component_construction_wall_ns = 0;
    uint64_t diameter_check_wall_ns = 0;
    // Time spent in online Dijkstra calls whose results depend only on the fixed
    // detector graph distances. Under a precomputed distance-lookup model, this
    // can be charged to per-code-distance preprocessing instead of per-shot ideal
    // parallel time.
    uint64_t precomputable_distance_lookup_wall_ns = 0;
    uint64_t influence_region_wall_ns = 0;
    uint64_t parent_assignment_wall_ns = 0;
    uint64_t growing_only_used_active_pair_table = 0;
    uint64_t growing_only_used_sparse_frontier = 0;
    uint64_t growing_only_used_component_lookup = 0;
    uint64_t growing_only_forced_max_level_shortcut_count = 0;
    uint64_t growing_only_processed_frontier_events = 0;
    uint64_t growing_only_generated_collision_events = 0;
    uint64_t growing_only_processed_collision_events = 0;
    uint64_t growing_only_exact_internal_pair_checks = 0;
    uint64_t growing_only_exact_external_pair_checks = 0;
    uint64_t growing_only_exact_diameter_checks = 0;
    uint64_t growing_only_max_component_active_size = 0;
    uint64_t growing_only_max_residual_active_size = 0;
};

struct ProcessingClusterRadiusNeighborList {
    cumulative_time_int buffer_bound = 0;
    std::vector<uint64_t> offsets;
    std::vector<uint32_t> neighbors;
};

struct ProcessingClusterGraphCache {
    size_t num_nodes = 0;
    bool has_all_pairs_distances = false;
    std::vector<uint32_t> packed_all_pairs_interior_distances;
    std::vector<uint64_t> packed_all_pairs_row_offsets;
    std::vector<cumulative_time_int> nearest_boundary_match_distance_by_vertex;
    std::vector<size_t> interior_component_id_by_vertex;
    std::vector<cumulative_time_int> interior_component_diameter_upper_bounds;
    std::vector<cumulative_time_int> boundary_edge_weight_by_vertex;
    std::vector<size_t> boundary_endpoint_vertices;
    cumulative_time_int max_finite_interior_distance = 0;
    cumulative_time_int interior_diameter_upper_bound = std::numeric_limits<cumulative_time_int>::max();
    cumulative_time_int max_boundary_edge_weight = 0;
    bool all_interior_pairs_finite = true;

    // Graph-only radius-neighbor lookup for fixed clustering parameters.
    // Levels are one-based; index 0 is intentionally unused.
    bool has_radius_neighbor_precompute = false;
    double radius_neighbor_beta = 0;
    double radius_neighbor_gamma = 0;
    double radius_neighbor_lambda = 0;
    size_t radius_neighbor_max_level = 0;
    std::vector<ProcessingClusterBounds> radius_neighbor_bounds_by_level;
    std::vector<ProcessingClusterRadiusNeighborList> radius_neighbors_by_level;
};


bool processing_cluster_bounds_phi_buffer_certificate_ok(
    const ProcessingClusterBounds& bounds);

ProcessingClusterBounds processing_cluster_bounds_for_level(
    size_t level, const ProcessingClusterConfig& config);

size_t processing_cluster_effective_max_level(const ProcessingClusterConfig& config);

std::vector<double> processing_cluster_stopping_lemma_phi_values(
    const ProcessingClusterConfig& config);

void validate_processing_cluster_stopping_lemma_bounds(
    const ProcessingClusterConfig& config);

std::vector<ProcessingClusterBounds> make_processing_cluster_stopping_lemma_bounds(
    const std::vector<cumulative_time_int>& diameter_bounds,
    cumulative_time_int min_buffer_slack = 1);

std::vector<ProcessingClusterBounds> make_processing_cluster_linear_stopping_lemma_bounds(
    cumulative_time_int diameter_unit,
    size_t num_levels,
    cumulative_time_int min_buffer_slack = 1);

// Converts an additive candidate diameter column into explicit clustering
// bounds.  The caller supplies candidates such as L_new, L_new + L_origin,
// L_new + 2 L_origin, ...; this function never changes the candidate column itself.
// When enforce_legacy_gap_filter is true, the selected subsequence skips ahead
// until d_candidate >= 3 d_prev + 4 b_prev before accepting the next level.
// When false, every usable additive candidate is considered consecutively.
// Buffers are chosen from the stopping-lemma phi recursion.  By default,
// min_buffer_to_diameter_ratio is also applied as a lower bound at every level.
// When min_buffer_ratio_first_level_only is true, that ratio is used only to
// choose b_1; levels k>1 are then controlled by the phi certificate, the
// minimal b_k>d_k condition, and monotonicity.
//
// Experimental phi-floor budget mode: when enforce_phi_floor_budget is true,
// only levels with phi_k >= phi_floor are selected, and b_k is additionally
// raised so the level's future phi-consumption T_k is at most
// phi_budget_per_level.  Combined with enforce_legacy_gap_filter, this keeps
// both phi and d_{k+1} >= 3 d_k + 4 b_k conditions explicit.
std::vector<ProcessingClusterBounds> make_processing_cluster_disjoint_stopping_lemma_bounds_from_candidates(
    const std::vector<cumulative_time_int>& candidate_diameter_bounds,
    size_t max_levels,
    cumulative_time_int min_buffer_slack = 1,
    double min_buffer_to_diameter_ratio = 1.0,
    bool enforce_legacy_gap_filter = false,
    bool min_buffer_ratio_first_level_only = false,
    bool enforce_phi_floor_budget = false,
    double phi_floor = 0.5,
    double phi_budget_per_level = 0.24);

// Direct-min-d variant of the disjoint stopping-lemma policy.  It does not
// build or scan an additive candidate column.  Instead it uses first_diameter_bound
// for level 1, then chooses the exact smallest integer next diameter satisfying
// d_{k+1} >= 3 d_k + 4 b_k for every later level.  Buffer and phi handling are
// otherwise the same as the phi-certified candidate policy above.
std::vector<ProcessingClusterBounds> make_processing_cluster_disjoint_stopping_lemma_bounds_direct_min_d(
    cumulative_time_int first_diameter_bound,
    size_t max_levels,
    cumulative_time_int min_buffer_slack = 1,
    double min_buffer_to_diameter_ratio = 1.0,
    bool min_buffer_ratio_first_level_only = false,
    bool enforce_phi_floor_budget = false,
    double phi_floor = 0.5,
    double phi_budget_per_level = 0.24);

// Target-level final-d optimizer.  It fixes the requested number of levels K,
// uses d_1=first_diameter_bound and d_{i+1}=3d_i+4b_i, then searches for a
// phi-consumption budget allocation that keeps phi_i >= phi_floor and minimizes
// the final diameter d_K.  This avoids additive candidate columns while also
// avoiding the fixed-per-level budget of the direct-min-d mode.
std::vector<ProcessingClusterBounds> make_processing_cluster_disjoint_stopping_lemma_bounds_optimize_final_d(
    cumulative_time_int first_diameter_bound,
    size_t target_levels,
    cumulative_time_int min_buffer_slack = 1,
    double min_buffer_to_diameter_ratio = 1.0,
    bool min_buffer_ratio_first_level_only = false,
    double phi_floor = 0.5,
    double phi_budget_total_fraction = 0.98,
    size_t search_rounds = 8);


// Phi-sequence direct-min-d policy.  It keeps the low-level buffers as small
// as possible while preserving a user-chosen decreasing lower envelope for
// the future stopping-lemma phi values:
//     phi_1 = 1,
//     target_phi_k = phi_floor + (phi2 - phi_floor) * q^(k-2), k>=2.
// Each b_k is chosen minimally so the actual next phi is at least
// target_phi_{k+1}; each next diameter is d_{k+1}=3d_k+4b_k.
std::vector<ProcessingClusterBounds> make_processing_cluster_disjoint_stopping_lemma_bounds_phi_sequence(
    cumulative_time_int first_diameter_bound,
    size_t max_levels,
    cumulative_time_int min_buffer_slack = 1,
    double min_buffer_to_diameter_ratio = 1.0,
    bool min_buffer_ratio_first_level_only = false,
    double phi2 = 0.0625,
    double phi_floor = 0.02,
    double phi_decay_q = 0.5,
    long double diameter_growth_additive = 0);

// Algorithm-style parameter schedule:
//     phi_1 = 1,
//     target_phi_{k+1} = phi_floor + (1-phi_floor) q^k,
//     d_{k+1} = 3 d_k + 4 b_k + 2 w_max.
// The runtime implementation uses d_1=w_max+1 and the integer-safe
// stopping-lemma inequalities
// 2(d_k+1)/(b_k-1) < phi_k and chooses the smallest certified integer buffer.
std::vector<ProcessingClusterBounds> make_processing_cluster_parameter_schedule_bounds(
    cumulative_time_int first_diameter_bound,
    size_t max_levels,
    double phi_floor = 0.02,
    double phi_decay_q = 0.5);

std::vector<size_t> make_processing_cluster_even_detector_count_caps(size_t num_levels);

ProcessingClusterGraphCache build_processing_cluster_graph_cache(const MatchingGraph& graph);

size_t infer_processing_cluster_precompute_max_level(
    const ProcessingClusterGraphCache& graph_cache,
    const ProcessingClusterConfig& config);

void ensure_processing_cluster_radius_neighbors_precomputed(
    const MatchingGraph& graph,
    ProcessingClusterGraphCache& graph_cache,
    const ProcessingClusterConfig& config);

ProcessingClusterHierarchy build_processing_cluster_hierarchy(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& active_detectors,
    const ProcessingClusterConfig& config = ProcessingClusterConfig(),
    const ProcessingClusterGraphCache* graph_cache = nullptr,
    ProcessingClusterProfilingStats* profiling_stats = nullptr);

ProcessingClusterHierarchy build_root_direct_processing_cluster_hierarchy(
    const MatchingGraph& graph,
    const std::vector<uint64_t>& active_detectors,
    const ProcessingClusterConfig& config = ProcessingClusterConfig(),
    const ProcessingClusterGraphCache* graph_cache = nullptr,
    ProcessingClusterProfilingStats* profiling_stats = nullptr);

}  // namespace pm

#endif  // PYMATCHING_PARALLEL_PROCESSING_CLUSTER_H
