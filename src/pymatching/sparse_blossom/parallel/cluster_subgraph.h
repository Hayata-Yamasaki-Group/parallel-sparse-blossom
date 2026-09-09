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

#ifndef PYMATCHING_PARALLEL_CLUSTER_SUBGRAPH_H
#define PYMATCHING_PARALLEL_CLUSTER_SUBGRAPH_H

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/user_graph.h"
#include "pymatching/sparse_blossom/parallel/processing_cluster.h"

namespace pm {

struct ClusterSubgraph {
    size_t cluster_id = std::numeric_limits<size_t>::max();
    UserGraph graph;
    // Fast path: prebuilt MWPM constructed directly from the already-discretized
    // global MatchingGraph. This avoids rebuilding/scaling a cluster UserGraph
    // during get_mwpm() in the lockstep scheduler. Only used for non-correlated
    // decoding with <=64 observables; other cases fall back to graph.get_mwpm().
    std::unique_ptr<Mwpm> cached_mwpm;
    // Optional non-owning MWPM for persistent full-graph worker states. When set,
    // ClusterSubgraph does not own or destroy the detector graph/MWPM state.
    // In the clean-pool mode this pointer is acquired from a ready queue before
    // the shot, used once, and then released to a background reset queue.
    Mwpm* persistent_mwpm = nullptr;
    // True when persistent_mwpm was acquired from the ready queue already reset.
    // The online setup should then skip an immediate reset; release will enqueue
    // the state for reset after extraction/export.
    bool persistent_mwpm_acquired_clean = false;
    // True for full-graph worker subgraphs whose persistent MWPM state should
    // be acquired lazily immediately before the cluster is initialized/executed.
    // This avoids holding one full-graph MWPM slot for every cluster during
    // per-shot subgraph construction.
    bool lazy_persistent_full_graph_worker_mwpm = false;
    // Optional non-owning MWPM for root full-graph identity clusters.  In the
    // canonical root=1 path this lets the root cluster reuse UserGraph's global
    // MWPM storage instead of cloning the whole MatchingGraph for every shot.
    // The scheduler resets this state after root extraction before the next shot.
    Mwpm* shared_mwpm = nullptr;
    std::vector<size_t> local_to_global_node_ids;
    std::vector<size_t> global_to_local_node_ids;
    std::vector<bool> local_boundary_node_mask;
    std::vector<uint64_t> global_active_detectors;
    std::vector<uint64_t> local_active_detectors;
    bool dropped_implied_weight_rules = false;
    // When true, the cluster subgraph is exactly the full global graph. The
    // level-batched scheduler may decode it using the already-built global
    // MWPM object instead of cloning the entire graph online. This preserves
    // the cluster/level structure; it only specializes the representation of a
    // full-graph influence subgraph.
    bool use_shared_full_graph_mwpm = false;
    // True when local detector ids are exactly the global detector ids.
    // Full-graph worker stopped-state imports use this to bypass all local/global remapping.
    bool local_node_ids_are_global_ids = false;
};

inline bool cluster_subgraph_local_node_is_boundary(const ClusterSubgraph& cluster_subgraph, size_t local_node) {
    return local_node < cluster_subgraph.local_boundary_node_mask.size() &&
        cluster_subgraph.local_boundary_node_mask[local_node];
}

ClusterSubgraph build_cluster_subgraph(UserGraph& graph, const ProcessingCluster& cluster);

// Builds a worker-local MWPM state by cloning the full detector MatchingGraph,
// while keeping cluster.active_detectors as the only input syndrome for that
// worker.  Unlike build_cluster_subgraph, this does not cut out an influence
// subgraph and therefore uses identity local/global detector ids.  The parent /
// child import machinery can still remap stopped states through these identity
// maps, preserving the cluster hierarchy without constructing per-cluster
// detector subgraphs.
ClusterSubgraph build_full_graph_worker_cluster_subgraph(
    UserGraph& graph,
    const ProcessingCluster& cluster,
    bool allow_shared_global_mwpm = true,
    bool lazy_acquire_persistent_mwpm = false);

// Pre-distributes full detector graph/MWPM states into persistent worker slots
// before online per-shot decoding starts. Later full-graph worker ClusterSubgraph
// objects acquire clean states non-owningly and release them for background reset.
// The number of states defaults to 4x worker_slots and can be overridden with
// PYMATCHING_PARENTLESS_STATE_POOL_SLOTS.
void prewarm_persistent_full_graph_worker_states(UserGraph& graph, size_t worker_slots);

// Waits until at least min_ready clean full-graph worker states are available.
// Used by the benchmark to separate warmup/pre-distribution from measured run time.
void wait_for_persistent_full_graph_worker_states_ready(UserGraph& graph, size_t min_ready);

// Returns a persistent state after use. In clean-pool mode this schedules reset
// on the background reset threads and immediately makes the caller stop using it.
void release_persistent_full_graph_worker_mwpm(Mwpm* mwpm);

// Acquires a clean persistent full-graph worker state for a lazy full-graph
// worker subgraph. No-op for non-lazy subgraphs and for already-acquired states.
void acquire_persistent_full_graph_worker_mwpm_for_subgraph(
    UserGraph& graph,
    ClusterSubgraph& cluster_subgraph);

Mwpm& cluster_subgraph_mwpm(ClusterSubgraph& cluster_subgraph, bool edge_correlations = false);

std::vector<uint64_t> map_global_detection_events_to_cluster_subgraph(
    const ClusterSubgraph& cluster_subgraph, const std::vector<uint64_t>& global_detection_events);

ExtendedMatchingResult decode_cluster_subgraph(
    ClusterSubgraph& cluster_subgraph,
    const std::vector<uint64_t>& global_detection_events,
    bool edge_correlations = false);

ExtendedMatchingResult decode_cluster_subgraph(
    ClusterSubgraph& cluster_subgraph, bool edge_correlations = false);

}  // namespace pm

#endif  // PYMATCHING_PARALLEL_CLUSTER_SUBGRAPH_H
