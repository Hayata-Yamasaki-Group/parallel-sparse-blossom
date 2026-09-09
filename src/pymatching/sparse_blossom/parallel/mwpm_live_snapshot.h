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

#ifndef PYMATCHING_PARALLEL_MWPM_LIVE_SNAPSHOT_H
#define PYMATCHING_PARALLEL_MWPM_LIVE_SNAPSHOT_H

#include <cstddef>
#include <vector>

#include "pymatching/sparse_blossom/flooder_matcher_interop/varying.h"
#include "pymatching/sparse_blossom/matcher/mwpm.h"
#include "pymatching/sparse_blossom/tracker/flood_check_event.h"

namespace pm {

constexpr size_t NO_LIVE_SNAPSHOT_NODE = SIZE_MAX;
constexpr size_t NO_LIVE_SNAPSHOT_REGION = SIZE_MAX;
constexpr size_t NO_LIVE_SNAPSHOT_ALT_TREE_NODE = SIZE_MAX;

struct LiveSnapshotQueuedEventTracker {
    cyclic_time_int desired_time{0};
    cyclic_time_int queued_time{0};
    bool has_desired_time = false;
    bool has_queued_time = false;
};

struct LiveSnapshotCompressedEdge {
    size_t loc_from = NO_LIVE_SNAPSHOT_NODE;
    size_t loc_to = NO_LIVE_SNAPSHOT_NODE;
    obs_int obs_mask = 0;
};

struct LiveSnapshotMatch {
    size_t region_id = NO_LIVE_SNAPSHOT_REGION;
    LiveSnapshotCompressedEdge edge;
};

struct LiveSnapshotRegionEdge {
    size_t region_id = NO_LIVE_SNAPSHOT_REGION;
    LiveSnapshotCompressedEdge edge;
};

struct DetectorNodeLiveSnapshot {
    size_t node_id = NO_LIVE_SNAPSHOT_NODE;
    size_t region_id = NO_LIVE_SNAPSHOT_REGION;
    size_t top_region_id = NO_LIVE_SNAPSHOT_REGION;
    size_t reached_from_source_id = NO_LIVE_SNAPSHOT_NODE;
    obs_int observables_crossed_from_source = 0;
    cumulative_time_int radius_of_arrival = 0;
    int32_t wrapped_radius_cached = 0;
    LiveSnapshotQueuedEventTracker node_event_tracker;
};

struct GraphFillRegionLiveSnapshot {
    size_t region_id = NO_LIVE_SNAPSHOT_REGION;
    size_t blossom_parent_id = NO_LIVE_SNAPSHOT_REGION;
    size_t blossom_parent_top_id = NO_LIVE_SNAPSHOT_REGION;
    size_t alt_tree_node_id = NO_LIVE_SNAPSHOT_ALT_TREE_NODE;
    VaryingCT radius;
    LiveSnapshotQueuedEventTracker shrink_event_tracker;
    LiveSnapshotMatch match;
    std::vector<LiveSnapshotRegionEdge> blossom_children;
    std::vector<size_t> shell_area_node_ids;
};

struct AltTreeEdgeLiveSnapshot {
    size_t alt_tree_node_id = NO_LIVE_SNAPSHOT_ALT_TREE_NODE;
    LiveSnapshotCompressedEdge edge;
};

struct AltTreeNodeLiveSnapshot {
    size_t alt_tree_node_id = NO_LIVE_SNAPSHOT_ALT_TREE_NODE;
    size_t inner_region_id = NO_LIVE_SNAPSHOT_REGION;
    size_t outer_region_id = NO_LIVE_SNAPSHOT_REGION;
    LiveSnapshotCompressedEdge inner_to_outer_edge;
    AltTreeEdgeLiveSnapshot parent;
    std::vector<AltTreeEdgeLiveSnapshot> children;
    bool visited = false;
};

struct FloodCheckEventLiveSnapshot {
    FloodCheckEventType type = NO_FLOOD_CHECK_EVENT;
    size_t target_id = SIZE_MAX;
    cyclic_time_int time{0};
};

struct FlooderQueueLiveSnapshot {
    cumulative_time_int cur_time = 0;
    std::vector<FloodCheckEventLiveSnapshot> events;
};

struct MwpmLiveSnapshot {
    size_t num_nodes = 0;
    size_t num_regions = 0;
    size_t num_alt_tree_nodes = 0;
    bool requires_search_graph = false;
    std::vector<DetectorNodeLiveSnapshot> detector_states;
    std::vector<GraphFillRegionLiveSnapshot> region_states;
    std::vector<AltTreeNodeLiveSnapshot> alt_tree_states;
    FlooderQueueLiveSnapshot flooder_queue;
    std::vector<CompressedEdge> match_edges;
    std::vector<uint64_t> negative_weight_detection_events;
    std::vector<size_t> negative_weight_observables;
    obs_int negative_weight_obs_mask = 0;
    total_weight_int negative_weight_sum = 0;
};

MwpmLiveSnapshot export_mwpm_live_snapshot(const Mwpm& mwpm);
void restore_mwpm_from_live_snapshot(Mwpm& target, const MwpmLiveSnapshot& snapshot);

}  // namespace pm

#endif  // PYMATCHING_PARALLEL_MWPM_LIVE_SNAPSHOT_H
