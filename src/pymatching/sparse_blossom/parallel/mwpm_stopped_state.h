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

#ifndef PYMATCHING_PARALLEL_MWPM_STOPPED_STATE_H
#define PYMATCHING_PARALLEL_MWPM_STOPPED_STATE_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include "pymatching/sparse_blossom/parallel/cluster_subgraph.h"
#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/flooder_matcher_interop/varying.h"

namespace pm {

constexpr size_t NO_STOPPED_NODE = SIZE_MAX;
constexpr size_t NO_STOPPED_REGION = SIZE_MAX;
constexpr size_t NO_STOPPED_ALT_TREE_NODE = SIZE_MAX;

struct StoppedCompressedEdge {
    size_t loc_from = NO_STOPPED_NODE;
    size_t loc_to = NO_STOPPED_NODE;
    obs_int obs_mask = 0;
};

struct StoppedMatch {
    size_t region_id = NO_STOPPED_REGION;
    StoppedCompressedEdge edge;
};

struct StoppedRegionEdge {
    size_t region_id = NO_STOPPED_REGION;
    StoppedCompressedEdge edge;
};

struct DetectorNodeStoppedState {
    size_t node_id = NO_STOPPED_NODE;
    size_t region_id = NO_STOPPED_REGION;
    size_t top_region_id = NO_STOPPED_REGION;
    size_t reached_from_source_id = NO_STOPPED_NODE;
    obs_int observables_crossed_from_source = 0;
    cumulative_time_int radius_of_arrival = 0;
    int32_t wrapped_radius_cached = 0;
};

struct GraphFillRegionStoppedState {
    size_t region_id = NO_STOPPED_REGION;
    size_t blossom_parent_id = NO_STOPPED_REGION;
    size_t blossom_parent_top_id = NO_STOPPED_REGION;
    size_t alt_tree_node_id = NO_STOPPED_ALT_TREE_NODE;
    VaryingCT radius;
    // Only exported for top-level stopped regions. Blossom descendants inherit their
    // local structure from blossom_children instead of a separately meaningful match.
    StoppedMatch match;
    std::vector<StoppedRegionEdge> blossom_children;
    std::vector<size_t> shell_area_node_ids;
};

struct StoppedStateImportOptions {
    bool reschedule_after_import = true;
    std::vector<size_t>* imported_target_local_node_ids = nullptr;

    // Experimental parentless/broadcast import mode.  A stopped state produced
    // by a lower-level cluster may be offered to an upper-level cluster that is
    // not its precomputed parent.  In that case parts of the stopped state can
    // lie outside the target cluster's buffer/influence subgraph, or can already
    // be present due to another broadcast path.  When this flag is set, import
    // first filters the source state down to the portion that can be safely
    // remapped into the target and does not overwrite already-owned detector
    // nodes.  If nothing remains, the import is a no-op.
    bool skip_unmappable_or_conflicting_state = false;

    // Experimental parentless/broadcast import mode for full-graph workers.
    // The same lower-level stopped state may be offered directly to an upper
    // cluster and then later arrive again through an intermediate level.  The
    // theory for the active-detector-limited export path says these duplicate
    // claims should be equivalent.  When this flag is set, detector ownership in
    // the target is overwritten by the imported state instead of rejecting or
    // filtering the import.
    bool allow_detector_ownership_overwrite = false;
};

struct MwpmStoppedState {
    cumulative_time_int algorithmic_time = 0;
    size_t num_nodes = 0;
    size_t num_regions = 0;
    size_t num_observables = 0;
    std::vector<DetectorNodeStoppedState> detector_states;
    std::vector<GraphFillRegionStoppedState> region_states;
    std::vector<uint64_t> negative_weight_detection_events;
    std::vector<size_t> negative_weight_observables;
    obs_int negative_weight_obs_mask = 0;
    total_weight_int negative_weight_sum = 0;
};

MwpmStoppedState export_mwpm_stopped_state(const Mwpm& mwpm);

// Experimental parentless export mode: keep only detector ownership whose
// flooder source is one of the cluster active detectors. The detector ids are
// local ids in the exporting cluster MWPM graph. This lets full-graph worker
// subgraphs broadcast/import stopped states without relying on a precomputed
// parent influence subgraph to crop the state.
MwpmStoppedState export_mwpm_stopped_state_for_active_detectors(
    const Mwpm& mwpm, const std::vector<uint64_t>& active_detector_local_node_ids);

MwpmStoppedState decode_detection_events_to_stopped_state(
    Mwpm& mwpm, const std::vector<uint64_t>& detection_events, bool edge_correlations = false);

std::vector<uint64_t> remap_source_detection_events_from_stopped_state(
    const ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    const MwpmStoppedState& source_state);

bool stopped_state_import_would_overwrite_detector_ownership(
    const ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    const MwpmStoppedState& source_state);

void import_stopped_state_into_cluster_subgraph(
    ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    const MwpmStoppedState& source_state,
    bool edge_correlations = false);

void import_stopped_state_into_cluster_subgraph(
    ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    const MwpmStoppedState& source_state,
    const StoppedStateImportOptions& options,
    bool edge_correlations = false);

void reschedule_imported_detector_nodes(
    ClusterSubgraph& target_subgraph,
    const std::vector<size_t>& target_local_node_ids,
    bool edge_correlations = false);

}  // namespace pm

#endif  // PYMATCHING_PARALLEL_MWPM_STOPPED_STATE_H
