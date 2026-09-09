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

#include "pymatching/sparse_blossom/parallel/mwpm_live_snapshot.h"

#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace pm {
namespace {

template <typename T>
std::vector<T*> collect_live_arena_objects(const Arena<T>& arena) {
    std::unordered_set<T*> available(arena.available.begin(), arena.available.end());
    std::vector<T*> live;
    live.reserve(arena.allocated.size() - arena.available.size());
    for (auto* ptr : arena.allocated) {
        if (!available.contains(ptr)) {
            live.push_back(ptr);
        }
    }
    return live;
}

LiveSnapshotQueuedEventTracker export_tracker(const QueuedEventTracker& tracker) {
    return {
        tracker.desired_time,
        tracker.queued_time,
        tracker.has_desired_time,
        tracker.has_queued_time,
    };
}

void restore_tracker(QueuedEventTracker& tracker, const LiveSnapshotQueuedEventTracker& snapshot) {
    tracker.desired_time = snapshot.desired_time;
    tracker.queued_time = snapshot.queued_time;
    tracker.has_desired_time = snapshot.has_desired_time;
    tracker.has_queued_time = snapshot.has_queued_time;
}

size_t node_id(const MatchingGraph& graph, const DetectorNode* node) {
    if (node == nullptr) {
        return NO_LIVE_SNAPSHOT_NODE;
    }
    auto id = static_cast<size_t>(node - graph.nodes.data());
    if (id >= graph.nodes.size()) {
        throw std::invalid_argument("Detector node pointer is outside the matching graph.");
    }
    return id;
}

LiveSnapshotCompressedEdge export_edge(const MatchingGraph& graph, const CompressedEdge& edge) {
    return {
        node_id(graph, edge.loc_from),
        node_id(graph, edge.loc_to),
        edge.obs_mask,
    };
}

CompressedEdge restore_edge(const MatchingGraph& graph, const LiveSnapshotCompressedEdge& edge) {
    auto* loc_from = edge.loc_from == NO_LIVE_SNAPSHOT_NODE ? nullptr : const_cast<DetectorNode*>(&graph.nodes[edge.loc_from]);
    auto* loc_to = edge.loc_to == NO_LIVE_SNAPSHOT_NODE ? nullptr : const_cast<DetectorNode*>(&graph.nodes[edge.loc_to]);
    return {loc_from, loc_to, edge.obs_mask};
}

size_t region_id(
    const GraphFillRegion* region,
    const std::unordered_map<const GraphFillRegion*, size_t>& region_ids) {
    if (region == nullptr) {
        return NO_LIVE_SNAPSHOT_REGION;
    }
    auto iter = region_ids.find(region);
    if (iter == region_ids.end()) {
        throw std::invalid_argument("GraphFillRegion pointer is not part of the live snapshot region arena.");
    }
    return iter->second;
}

size_t alt_tree_id(
    const AltTreeNode* node,
    const std::unordered_map<const AltTreeNode*, size_t>& node_ids) {
    if (node == nullptr) {
        return NO_LIVE_SNAPSHOT_ALT_TREE_NODE;
    }
    auto iter = node_ids.find(node);
    if (iter == node_ids.end()) {
        throw std::invalid_argument("AltTreeNode pointer is not part of the live snapshot node arena.");
    }
    return iter->second;
}

FloodCheckEventLiveSnapshot export_flood_check_event(
    const MatchingGraph& graph,
    const FloodCheckEvent& event,
    const std::unordered_map<const GraphFillRegion*, size_t>& region_ids) {
    switch (event.tentative_event_type) {
        case NO_FLOOD_CHECK_EVENT:
            return {NO_FLOOD_CHECK_EVENT, SIZE_MAX, event.time};
        case LOOK_AT_NODE:
            return {LOOK_AT_NODE, node_id(graph, event.data_look_at_node), event.time};
        case LOOK_AT_SHRINKING_REGION:
            return {LOOK_AT_SHRINKING_REGION, region_id(event.data_look_at_shrinking_region, region_ids), event.time};
        case LOOK_AT_SEARCH_NODE:
            throw std::invalid_argument(
                "Live MWPM snapshots do not yet support active SearchFlooder events.");
    }
    throw std::invalid_argument("Unrecognized FloodCheckEventType.");
}

FloodCheckEvent restore_flood_check_event(
    Mwpm& mwpm,
    const FloodCheckEventLiveSnapshot& event,
    const std::vector<GraphFillRegion*>& regions) {
    switch (event.type) {
        case NO_FLOOD_CHECK_EVENT:
            return FloodCheckEvent(event.time);
        case LOOK_AT_NODE:
            return FloodCheckEvent(&mwpm.flooder.graph.nodes[event.target_id], event.time);
        case LOOK_AT_SHRINKING_REGION:
            return FloodCheckEvent(regions[event.target_id], event.time);
        case LOOK_AT_SEARCH_NODE:
            throw std::invalid_argument(
                "Live MWPM snapshots do not yet support restoring active SearchFlooder events.");
    }
    throw std::invalid_argument("Unrecognized FloodCheckEventType.");
}

template <typename T>
void reset_arena_to_empty(Arena<T>& arena) {
    arena.~Arena<T>();
    new (&arena) Arena<T>();
}

void clear_mwpm_live_state(Mwpm& mwpm) {
    for (auto& node : mwpm.flooder.graph.nodes) {
        node.reset();
    }
    mwpm.flooder.queue.reset();
    mwpm.flooder.match_edges.clear();
    mwpm.flooder.negative_weight_detection_events.clear();
    mwpm.flooder.negative_weight_observables.clear();
    mwpm.flooder.negative_weight_obs_mask = 0;
    mwpm.flooder.negative_weight_sum = 0;
    reset_arena_to_empty(mwpm.flooder.region_arena);
    reset_arena_to_empty(mwpm.node_arena);
    mwpm.search_flooder.reset();
    mwpm.search_flooder.target_type = NO_TARGET;
}

}  // namespace

MwpmLiveSnapshot export_mwpm_live_snapshot(const Mwpm& mwpm) {
    if (!mwpm.search_flooder.queue.empty() || !mwpm.search_flooder.reached_nodes.empty() ||
        mwpm.search_flooder.target_type != NO_TARGET) {
        throw std::invalid_argument(
            "Live MWPM snapshots do not yet support non-idle SearchFlooder state.");
    }

    MwpmLiveSnapshot snapshot;
    snapshot.num_nodes = mwpm.flooder.graph.nodes.size();
    snapshot.requires_search_graph = !mwpm.search_flooder.graph.nodes.empty();
    snapshot.match_edges = mwpm.flooder.match_edges;
    snapshot.negative_weight_detection_events = mwpm.flooder.negative_weight_detection_events;
    snapshot.negative_weight_observables = mwpm.flooder.negative_weight_observables;
    snapshot.negative_weight_obs_mask = mwpm.flooder.negative_weight_obs_mask;
    snapshot.negative_weight_sum = mwpm.flooder.negative_weight_sum;

    auto live_regions = collect_live_arena_objects(mwpm.flooder.region_arena);
    auto live_alt_tree_nodes = collect_live_arena_objects(mwpm.node_arena);
    snapshot.num_regions = live_regions.size();
    snapshot.num_alt_tree_nodes = live_alt_tree_nodes.size();

    std::unordered_map<const GraphFillRegion*, size_t> region_ids;
    for (size_t k = 0; k < live_regions.size(); k++) {
        region_ids.emplace(live_regions[k], k);
    }
    std::unordered_map<const AltTreeNode*, size_t> alt_tree_ids;
    for (size_t k = 0; k < live_alt_tree_nodes.size(); k++) {
        alt_tree_ids.emplace(live_alt_tree_nodes[k], k);
    }

    snapshot.detector_states.reserve(mwpm.flooder.graph.nodes.size());
    for (size_t k = 0; k < mwpm.flooder.graph.nodes.size(); k++) {
        const auto& node = mwpm.flooder.graph.nodes[k];
        snapshot.detector_states.push_back({
            k,
            region_id(node.region_that_arrived, region_ids),
            region_id(node.region_that_arrived_top, region_ids),
            node_id(mwpm.flooder.graph, node.reached_from_source),
            node.observables_crossed_from_source,
            node.radius_of_arrival,
            node.wrapped_radius_cached,
            export_tracker(node.node_event_tracker),
        });
    }

    snapshot.region_states.reserve(live_regions.size());
    for (size_t k = 0; k < live_regions.size(); k++) {
        const auto* region = live_regions[k];
        GraphFillRegionLiveSnapshot region_state;
        region_state.region_id = k;
        region_state.blossom_parent_id = region_id(region->blossom_parent, region_ids);
        region_state.blossom_parent_top_id = region_id(region->blossom_parent_top, region_ids);
        region_state.alt_tree_node_id = alt_tree_id(region->alt_tree_node, alt_tree_ids);
        region_state.radius = region->radius;
        region_state.shrink_event_tracker = export_tracker(region->shrink_event_tracker);
        if (region->blossom_parent == nullptr && region->alt_tree_node == nullptr) {
            region_state.match = {
                region_id(region->match.region, region_ids),
                export_edge(mwpm.flooder.graph, region->match.edge),
            };
        }
        region_state.blossom_children.reserve(region->blossom_children.size());
        for (const auto& child : region->blossom_children) {
            region_state.blossom_children.push_back({
                region_id(child.region, region_ids),
                export_edge(mwpm.flooder.graph, child.edge),
            });
        }
        region_state.shell_area_node_ids.reserve(region->shell_area.size());
        for (const auto* node : region->shell_area) {
            region_state.shell_area_node_ids.push_back(node_id(mwpm.flooder.graph, node));
        }
        snapshot.region_states.push_back(std::move(region_state));
    }

    snapshot.alt_tree_states.reserve(live_alt_tree_nodes.size());
    for (size_t k = 0; k < live_alt_tree_nodes.size(); k++) {
        const auto* node = live_alt_tree_nodes[k];
        AltTreeNodeLiveSnapshot alt_tree_state;
        alt_tree_state.alt_tree_node_id = k;
        alt_tree_state.inner_region_id = region_id(node->inner_region, region_ids);
        alt_tree_state.outer_region_id = region_id(node->outer_region, region_ids);
        alt_tree_state.inner_to_outer_edge = export_edge(mwpm.flooder.graph, node->inner_to_outer_edge);
        alt_tree_state.parent = {
            alt_tree_id(node->parent.alt_tree_node, alt_tree_ids),
            export_edge(mwpm.flooder.graph, node->parent.edge),
        };
        alt_tree_state.children.reserve(node->children.size());
        for (const auto& child : node->children) {
            alt_tree_state.children.push_back({
                alt_tree_id(child.alt_tree_node, alt_tree_ids),
                export_edge(mwpm.flooder.graph, child.edge),
            });
        }
        alt_tree_state.visited = node->visited;
        snapshot.alt_tree_states.push_back(std::move(alt_tree_state));
    }

    snapshot.flooder_queue.cur_time = mwpm.flooder.queue.cur_time;
    for (const auto& bucket : mwpm.flooder.queue.bit_buckets) {
        for (const auto& event : bucket) {
            snapshot.flooder_queue.events.push_back(
                export_flood_check_event(mwpm.flooder.graph, event, region_ids));
        }
    }

    return snapshot;
}

void restore_mwpm_from_live_snapshot(Mwpm& restored, const MwpmLiveSnapshot& snapshot) {
    if (restored.flooder.graph.nodes.size() != snapshot.num_nodes) {
        throw std::invalid_argument("Target MWPM graph is incompatible with live snapshot node count.");
    }
    bool target_has_search_graph = !restored.search_flooder.graph.nodes.empty();
    if (target_has_search_graph != snapshot.requires_search_graph) {
        throw std::invalid_argument("Target MWPM search-graph presence is incompatible with live snapshot.");
    }
    clear_mwpm_live_state(restored);

    restored.flooder.negative_weight_detection_events = snapshot.negative_weight_detection_events;
    restored.flooder.negative_weight_observables = snapshot.negative_weight_observables;
    restored.flooder.negative_weight_obs_mask = snapshot.negative_weight_obs_mask;
    restored.flooder.negative_weight_sum = snapshot.negative_weight_sum;
    restored.flooder.match_edges = snapshot.match_edges;

    std::vector<GraphFillRegion*> regions(snapshot.region_states.size(), nullptr);
    for (const auto& region_state : snapshot.region_states) {
        auto* region = restored.flooder.region_arena.alloc_default_constructed();
        region->radius = region_state.radius;
        restore_tracker(region->shrink_event_tracker, region_state.shrink_event_tracker);
        region->match.clear();
        region->blossom_children.clear();
        region->shell_area.clear();
        regions[region_state.region_id] = region;
    }

    std::vector<AltTreeNode*> alt_tree_nodes(snapshot.alt_tree_states.size(), nullptr);
    for (const auto& alt_tree_state : snapshot.alt_tree_states) {
        auto* node = restored.node_arena.alloc_default_constructed();
        alt_tree_nodes[alt_tree_state.alt_tree_node_id] = node;
    }

    for (const auto& region_state : snapshot.region_states) {
        auto* region = regions[region_state.region_id];
        region->blossom_parent = region_state.blossom_parent_id == NO_LIVE_SNAPSHOT_REGION
            ? nullptr
            : regions[region_state.blossom_parent_id];
        region->blossom_parent_top = region_state.blossom_parent_top_id == NO_LIVE_SNAPSHOT_REGION
            ? nullptr
            : regions[region_state.blossom_parent_top_id];
        region->alt_tree_node = region_state.alt_tree_node_id == NO_LIVE_SNAPSHOT_ALT_TREE_NODE
            ? nullptr
            : alt_tree_nodes[region_state.alt_tree_node_id];
        region->match = {
            region_state.match.region_id == NO_LIVE_SNAPSHOT_REGION ? nullptr : regions[region_state.match.region_id],
            restore_edge(restored.flooder.graph, region_state.match.edge),
        };
        for (const auto& child : region_state.blossom_children) {
            region->blossom_children.push_back({
                child.region_id == NO_LIVE_SNAPSHOT_REGION ? nullptr : regions[child.region_id],
                restore_edge(restored.flooder.graph, child.edge),
            });
        }
        for (auto node_id_value : region_state.shell_area_node_ids) {
            region->shell_area.push_back(&restored.flooder.graph.nodes[node_id_value]);
        }
    }

    for (const auto& alt_tree_state : snapshot.alt_tree_states) {
        auto* node = alt_tree_nodes[alt_tree_state.alt_tree_node_id];
        node->inner_region = alt_tree_state.inner_region_id == NO_LIVE_SNAPSHOT_REGION
            ? nullptr
            : regions[alt_tree_state.inner_region_id];
        node->outer_region = alt_tree_state.outer_region_id == NO_LIVE_SNAPSHOT_REGION
            ? nullptr
            : regions[alt_tree_state.outer_region_id];
        node->inner_to_outer_edge = restore_edge(restored.flooder.graph, alt_tree_state.inner_to_outer_edge);
        node->parent = {
            alt_tree_state.parent.alt_tree_node_id == NO_LIVE_SNAPSHOT_ALT_TREE_NODE
                ? nullptr
                : alt_tree_nodes[alt_tree_state.parent.alt_tree_node_id],
            restore_edge(restored.flooder.graph, alt_tree_state.parent.edge),
        };
        node->children.clear();
        for (const auto& child : alt_tree_state.children) {
            node->children.push_back({
                child.alt_tree_node_id == NO_LIVE_SNAPSHOT_ALT_TREE_NODE ? nullptr : alt_tree_nodes[child.alt_tree_node_id],
                restore_edge(restored.flooder.graph, child.edge),
            });
        }
        node->visited = alt_tree_state.visited;
    }

    for (const auto& detector_state : snapshot.detector_states) {
        auto& node = restored.flooder.graph.nodes[detector_state.node_id];
        node.region_that_arrived = detector_state.region_id == NO_LIVE_SNAPSHOT_REGION
            ? nullptr
            : regions[detector_state.region_id];
        node.region_that_arrived_top = detector_state.top_region_id == NO_LIVE_SNAPSHOT_REGION
            ? nullptr
            : regions[detector_state.top_region_id];
        node.reached_from_source = detector_state.reached_from_source_id == NO_LIVE_SNAPSHOT_NODE
            ? nullptr
            : &restored.flooder.graph.nodes[detector_state.reached_from_source_id];
        node.observables_crossed_from_source = detector_state.observables_crossed_from_source;
        node.radius_of_arrival = detector_state.radius_of_arrival;
        node.wrapped_radius_cached = detector_state.wrapped_radius_cached;
        restore_tracker(node.node_event_tracker, detector_state.node_event_tracker);
    }

    restored.flooder.queue.reset();
    restored.flooder.queue.cur_time = snapshot.flooder_queue.cur_time;
    for (const auto& event : snapshot.flooder_queue.events) {
        if (event.type == NO_FLOOD_CHECK_EVENT) {
            continue;
        }
        restored.flooder.queue.enqueue(restore_flood_check_event(restored, event, regions));
    }

}

}  // namespace pm
