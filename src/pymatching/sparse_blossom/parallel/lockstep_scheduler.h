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

#ifndef PYMATCHING_PARALLEL_LOCKSTEP_SCHEDULER_H
#define PYMATCHING_PARALLEL_LOCKSTEP_SCHEDULER_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/user_graph.h"
#include "pymatching/sparse_blossom/parallel/cluster_subgraph.h"
#include "pymatching/sparse_blossom/parallel/mwpm_live_snapshot.h"
#include "pymatching/sparse_blossom/parallel/mwpm_stopped_state.h"
#include "pymatching/sparse_blossom/parallel/processing_cluster.h"

namespace pm {

enum class LockstepClusterExecutionState : uint8_t {
    RUNNING,
    QUIESCENT_STOPPED,
    DRAINED_WITH_LIVE_ALT_TREE,
    IMPORTED_INTO_PARENT,
};

struct LockstepSchedulerConfig {
    ProcessingClusterConfig cluster_config;
    bool edge_correlations = false;
    bool allow_global_fallback = true;
    size_t num_workers = 1;
    // Experimental safe/debug mode: keep overlap checkpoints for recovering parent
    // state when a stopped child would overwrite detector ownership.
    //
    // The normal fast measurement path leaves this off.  If an overwrite conflict
    // is encountered with checkpoints disabled, the parent cluster is rebuilt from
    // algorithmic time 0 and already-imported siblings are replayed in stop-time
    // order before importing the new child.  The checkpoint code is intentionally
    // retained for experiments and can be re-enabled by callers.
    bool enable_overlap_checkpoints = false;
    // Experimental fast benchmark mode: if the hierarchy has exactly one root cluster, decode that
    // root subgraph directly using all input detection events that map into it instead of running
    // the full child/parent lockstep scheduler. This is intentionally opt-in because it bypasses
    // the paper-style state inheritance machinery; callers should compare predicted observables
    // against the global decoder when using it.
    bool enable_single_root_direct_decode_fast_path = false;
    // More general version of the benchmark shortcut: decode only the root / maximum-level
    // processing clusters directly and combine their observable flips. This removes lower-level
    // scheduler/inheritance overhead and measures the paper-style ideal unit of work more closely.
    // It is only safe to use after comparing against the global decoder for the target workload.
    bool enable_root_direct_decode_fast_path = false;
    // Even more aggressive verified benchmarking mode: after constructing root
    // coverage clusters, split the detection events by root influence region but
    // decode each split on the original global graph. This bypasses the
    // parent/child import machinery and is therefore not the default.
    bool enable_root_global_decode_fast_path = false;
    // No-cutout worker graph mode: preserve the processing-cluster hierarchy and
    // stopped-state parent/child import, but build every worker MWPM on a full
    // copy of the original detector graph instead of constructing a restricted
    // ClusterSubgraph/influence subgraph.  Local detector ids are identity-mapped
    // to global detector ids for import/export.
    bool enable_full_graph_worker_subgraphs = false;
    // Accounting model for enable_full_graph_worker_subgraphs: the full detector
    // graph is distributed to workers before the online decoding window starts.
    // Therefore worker-graph materialization/clone time is reported separately
    // and excluded from preprocessing, initial setup, and the critical path.
    bool exclude_full_graph_worker_distribution_from_runtime = true;
    // Paper-style barrier scheduler: process clusters level by level, importing
    // stopped lower-level states into parents before processing the next level.
    // This preserves the hierarchy/inheritance structure while avoiding the
    // expensive global event-time lockstep scheduler.
    bool enable_level_batched_scheduler = false;
    // Experimental import/export policy that does not use precomputed
    // parent-child edges as the sole import route.  When a lower-level cluster
    // stops, its stopped state is first filtered on export to ownership reached
    // from that cluster's active detectors, then offered to every higher-level
    // non-empty cluster.  This is intended for clustering paths that do not
    // materialize influence subgraphs; upper full-graph workers can import the
    // already-filtered state without knowing a parent-child relation.
    bool enable_parentless_broadcast_import = false;
    // Experimental scheduler integration for the fast growing-only clustering
    // path.  Instead of building influence regions and parent-child edges with
    // build_processing_cluster_hierarchy, construct the online hierarchy from
    // GrowingOnlyClusterer results.  This is intended to be used together with
    // full-graph workers and parentless broadcast import.
    bool enable_growing_only_processing_clusters = false;
    // Parentless-specific level pipeline.  For growing-only clusters there is no
    // parent/child dependency graph to maintain: execute one whole level, export
    // the stopped states from that level, broadcast them to upper non-empty
    // clusters with overwrite enabled, then continue.  This bypasses the old
    // lockstep scheduler's global-time scans, per-cluster readiness checks, and
    // parent lookup bookkeeping.
    bool enable_parentless_level_pipeline = false;
    // Legacy fast path for the degenerate case where every cluster is marked as a
    // root. It decodes all clusters directly and therefore cannot implement the
    // parentless broadcast critical-event accounting based on previous-level
    // cutoff times. Keep it opt-in and disabled for event-count comparisons.
    bool enable_parentless_all_root_direct_io_path = false;
    // Worker count used only for ideal multi-core accounting. The actual benchmark
    // may still run with num_workers=1 so that preprocessing remains on one core
    // and per-cluster timings are measured without thread-dispatch noise. A value
    // of 0 means use num_workers.
    size_t ideal_cluster_worker_count = 0;

    // Optional externally-built hierarchy. When supplied, the scheduler skips its
    // normal clustering stage and executes this hierarchy with the same
    // parentless import/export pipeline. This is used to compare different
    // clustering policies while keeping the sparse-blossom event executor
    // identical. The caller must keep the pointed-to hierarchy alive for the
    // duration of lockstep_hierarchical_decode().
    const ProcessingClusterHierarchy* prebuilt_processing_cluster_hierarchy = nullptr;
};

struct LockstepProfilingStats {
    uint64_t clustering_wall_ns = 0;
    uint64_t clustering_component_construction_wall_ns = 0;
    uint64_t clustering_diameter_check_wall_ns = 0;
    uint64_t clustering_precomputable_distance_lookup_wall_ns = 0;
    uint64_t clustering_influence_region_wall_ns = 0;
    uint64_t clustering_parent_assignment_wall_ns = 0;
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
    uint64_t subgraph_build_wall_ns = 0;
    uint64_t worker_graph_distribution_wall_ns = 0;
    uint64_t single_root_decode_wall_ns = 0;
    uint64_t root_extraction_wall_ns = 0;
    uint64_t initial_mwpm_setup_wall_ns = 0;
    uint64_t scheduler_wall_ns = 0;
    uint64_t child_import_wall_ns = 0;
    uint64_t overlap_checkpoint_capture_wall_ns = 0;
    uint64_t overlap_checkpoint_restore_wall_ns = 0;
    size_t cluster_count = 0;
    size_t root_cluster_count = 0;
    size_t max_level = 0;
    size_t max_clusters_in_level = 0;
    size_t sum_cluster_active_detectors = 0;
    size_t max_cluster_active_detectors = 0;
    size_t child_import_count = 0;
    size_t overlap_checkpoint_capture_count = 0;
    size_t overlap_checkpoint_restore_count = 0;
    size_t scheduler_batch_count = 0;
    size_t scheduler_cluster_step_count = 0;
    size_t import_candidate_count = 0;
    uint64_t cluster_step_wall_ns = 0;
    uint64_t max_cluster_step_wall_ns = 0;
    uint64_t max_level_cluster_step_wall_ns = 0;

    // Direct all-root worker breakdown. These timers split the direct worker
    // path into the stages that are hidden inside
    // decode_detection_events_for_up_to_64_observables(...):
    //   init:    active detector insertion into the worker MWPM state
    //   advance: sparse-blossom timeline processing after initialization
    //   extract: compact observable extraction from the quiescent state
    //   release: returning/resetting the persistent state slot
    // Sums are over all root workers. Max values are the actual worker critical
    // path for each stage.
    uint64_t direct_worker_init_wall_ns = 0;
    uint64_t direct_worker_advance_wall_ns = 0;
    uint64_t direct_worker_extract_wall_ns = 0;
    uint64_t direct_worker_release_wall_ns = 0;
    uint64_t direct_worker_total_wall_ns = 0;
    uint64_t max_direct_worker_init_wall_ns = 0;
    uint64_t max_direct_worker_advance_wall_ns = 0;
    uint64_t max_direct_worker_extract_wall_ns = 0;
    uint64_t max_direct_worker_release_wall_ns = 0;
    uint64_t max_direct_worker_total_wall_ns = 0;
    uint64_t direct_worker_active_detectors = 0;
    uint64_t max_direct_worker_active_detectors = 0;

    // Sparse-blossom MWPM event counts measured from scheduler-driven worker
    // states, after child stopped-state imports/reschedules have happened.
    // parallel_cluster_events is the serial total over all raw cluster advances.
    // In the parentless level pipeline, parallel_level_max_event_sum keeps its
    // historical field name but now stores the level-by-level ideal critical-path
    // event count:
    //   cumulative_events[level] = cumulative_events[level - 1]
    //       + max_cluster(events occurring after previous_level_stop_time).
    // Events at or before the previous level's maximum stopping time are treated
    // as overlappable with lower-level work and are not added. The same cutoff is
    // applied to algorithmic stopping times. For non-pipeline paths, the old
    // max-cluster fallback is retained.
    uint64_t parallel_cluster_events = 0;
    uint64_t parallel_level_max_event_sum = 0;
    uint64_t parallel_max_cluster_events = 0;
    uint64_t parallel_nonempty_event_levels = 0;
    cumulative_time_int parallel_level_critical_path_algorithmic_time = 0;

    // Parallel timing policy used by local_parallel_bench:
    //
    // We assume clusters from all levels can be executed concurrently and that
    // each worker already owns its detector graph/MWPM state before online
    // decoding starts. Therefore the comparison metric is the sparse-blossom
    // execution critical path of the slowest cluster over all levels, not the
    // serial scheduler wall time and not a sum over levels. It includes imports
    // and scheduler-driven execution work performed by that cluster:
    //
    //   event_step + child_import
    //     + checkpoint_capture + checkpoint_restore + mark_stop_or_export
    //
    // It intentionally excludes detector-graph distribution, cut-out subgraph
    // materialization, and initial MWPM setup. The excluded setup for the same
    // selected critical-path cluster is still reported separately for diagnosis.
    uint64_t max_level_cluster_lifecycle_wall_ns = 0;
    uint64_t max_level_cluster_initial_setup_wall_ns = 0;
    uint64_t max_level_cluster_child_import_wall_ns = 0;
    uint64_t max_level_cluster_checkpoint_capture_wall_ns = 0;
    uint64_t max_level_cluster_checkpoint_restore_wall_ns = 0;
    uint64_t max_level_cluster_mark_stop_wall_ns = 0;
    uint64_t max_level_cluster_id = 0;
    // Estimated sparse-blossom worker time under a persistent worker pool. This is
    // the greedy/LPT makespan over all cluster_step_wall_ns values at once, not a
    // sum over levels. With enough workers for all clusters, this reduces to the
    // slowest cluster over all levels.
    uint64_t ideal_worker_cluster_makespan_wall_ns = 0;
    // Simple shared-memory communication estimate for the ideal worker-pool model.
    // This includes level barriers, per-cluster dispatch/inbox bookkeeping, child
    // import messages, and root-result reduction.
    uint64_t ideal_worker_communication_estimate_ns = 0;
    size_t ideal_cluster_worker_count_used = 0;
};

struct LockstepOverlapSnapshotHistoryEntry {
    cumulative_time_int time = 0;
    MwpmLiveSnapshot snapshot;
    std::vector<uint64_t> extraction_detection_events;
};

struct LockstepDecodeResult {
    ProcessingClusterHierarchy hierarchy;
    std::vector<ClusterSubgraph> cluster_subgraphs;
    std::vector<LockstepClusterExecutionState> cluster_execution_states;
    std::vector<MwpmStoppedState> cluster_stopped_states;
    std::vector<bool> cluster_has_stopped_state;
    std::vector<bool> cluster_result_is_provisional;
    std::vector<bool> cluster_has_overlap_checkpoint;
    std::vector<MwpmLiveSnapshot> cluster_overlap_checkpoints;
    std::vector<cumulative_time_int> cluster_overlap_checkpoint_times;
    std::vector<std::vector<uint64_t>> cluster_overlap_checkpoint_extraction_detection_events;
    std::vector<bool> cluster_used_overlap_checkpoint_recovery;
    // Per-cluster live-state snapshot history.  When a child stops after its
    // parent has already advanced, restore the parent to the newest snapshot at
    // or before the child stop time, then replay already-imported siblings in
    // chronological order.  This avoids falling back to the time-0 snapshot.
    std::vector<std::vector<LockstepOverlapSnapshotHistoryEntry>> cluster_overlap_snapshot_histories;
    // Precomputed scheduler topology. direct_child_cluster_ids avoids scanning all
    // clusters whenever a parent imports/checkpoints children.
    std::vector<std::vector<size_t>> direct_child_cluster_ids;
    // For each parent cluster and parent-local detector node, the child clusters whose
    // boundary/non-boundary nodes would be affected if that detector were claimed.
    // This replaces per-event per-child watched-mask construction.
    std::vector<std::vector<std::vector<size_t>>> watcher_child_ids_by_parent_local_node;
    std::vector<cumulative_time_int> cluster_last_progress_time;
    std::vector<uint64_t> cluster_step_wall_ns;
    std::vector<uint64_t> cluster_initial_setup_wall_ns;
    std::vector<uint64_t> cluster_child_import_wall_ns;
    std::vector<uint64_t> cluster_checkpoint_capture_wall_ns;
    std::vector<uint64_t> cluster_checkpoint_restore_wall_ns;
    std::vector<uint64_t> cluster_mark_stop_wall_ns;
    // Raw MWPM notifications actually processed while advancing this cluster.
    std::vector<size_t> cluster_step_event_counts;
    // Level-by-level ideal critical-path accounting for the parentless pipeline.
    // For a cluster in level k>1, events at or before the previous level's
    // maximum algorithmic stop time can overlap lower-level execution.  The
    // overlap phase still costs the maximum of (the critical events carried from
    // lower levels) and (this cluster's events at/before the cutoff).  Events
    // strictly after the cutoff are then added for this cluster.  The level's
    // carried critical event count is the maximum cluster critical count in that
    // level.
    std::vector<size_t> cluster_critical_path_event_counts;
    std::vector<uint8_t> cluster_has_critical_path_event_count;
    // Per-cluster event-accounting details used to validate the parentless
    // level pipeline.  A cutoff exists for levels above the first completed
    // level and is the previous completed level's maximum stop time.
    std::vector<uint8_t> cluster_has_cutoff_algorithmic_time;
    std::vector<cumulative_time_int> cluster_cutoff_algorithmic_times;
    std::vector<size_t> cluster_prior_level_critical_event_counts;
    std::vector<size_t> cluster_overlap_event_bases;
    std::vector<size_t> cluster_events_at_or_before_cutoff;
    std::vector<size_t> cluster_events_after_cutoff;
    std::vector<cumulative_time_int> cluster_stop_algorithmic_times;
    std::vector<cumulative_time_int> cluster_critical_path_algorithmic_times;
    std::vector<std::vector<uint64_t>> cluster_extraction_detection_events;
    // Per-cluster marker table for cluster_extraction_detection_events. This avoids
    // sort/unique after every child import; merging a child becomes O(new sources).
    std::vector<std::vector<uint8_t>> cluster_extraction_detection_event_seen;
    // Root clusters are not imported into a parent.  For <=64 observables, keep their
    // final matching in the same compact obs_mask + weight form used by the global
    // sparse-blossom fast path.  This avoids expanding into ExtendedMatchingResult
    // during the root critical path.  The ExtendedMatchingResult path is retained for
    // experiments and for graphs with more than 64 observables.
    std::vector<MatchingResult> cluster_direct_compact_results;
    std::vector<bool> cluster_has_direct_compact_result;
    // Debug-only per-cluster result snapshots. These are populated for both
    // terminal/direct clusters and provisional stopped-state exports, but are
    // never used for root aggregation. They are intended for mistake logs.
    std::vector<MatchingResult> cluster_debug_compact_results;
    std::vector<bool> cluster_has_debug_compact_result;
    std::vector<ExtendedMatchingResult> cluster_direct_results;
    std::vector<bool> cluster_has_direct_result;
    ExtendedMatchingResult root_aggregate_result;
    LockstepProfilingStats profiling_stats;
    bool used_global_fallback = false;
    bool encountered_non_quiescent_cluster = false;
};

LockstepDecodeResult lockstep_hierarchical_decode(
    UserGraph& graph,
    const std::vector<uint64_t>& active_detectors,
    const LockstepSchedulerConfig& config = LockstepSchedulerConfig());

// Constructs and wakes the persistent worker pool used by the parentless
// level-pipeline actual parallel benchmark. The benchmark driver calls this
// before the online per-shot timer so thread creation is not charged to the
// measured decoding window. It is a no-op for worker_count <= 1.
void prewarm_parentless_level_worker_pool(size_t worker_count);

}  // namespace pm

#endif  // PYMATCHING_PARALLEL_LOCKSTEP_SCHEDULER_H
