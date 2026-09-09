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

#include "pymatching/sparse_blossom/parallel/mwpm_stopped_state.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>

#include "pymatching/sparse_blossom/flooder/detector_node.h"
#include "pymatching/sparse_blossom/flooder/graph_fill_region.h"

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

size_t stopped_node_id(const MatchingGraph& graph, const DetectorNode* node) {
    if (node == nullptr) {
        return NO_STOPPED_NODE;
    }
    auto node_id = static_cast<size_t>(node - graph.nodes.data());
    if (node_id >= graph.nodes.size()) {
        throw std::invalid_argument("DetectorNode pointer is outside Mwpm flooder graph.");
    }
    return node_id;
}

size_t stopped_region_id(
    const GraphFillRegion* region, const std::unordered_map<const GraphFillRegion*, size_t>& region_ids) {
    if (region == nullptr) {
        return NO_STOPPED_REGION;
    }
    auto iter = region_ids.find(region);
    if (iter == region_ids.end()) {
        throw std::invalid_argument("GraphFillRegion pointer is not part of the stopped-state region arena.");
    }
    return iter->second;
}

StoppedCompressedEdge export_stopped_compressed_edge(const MatchingGraph& graph, const CompressedEdge& edge) {
    return {
        stopped_node_id(graph, edge.loc_from),
        stopped_node_id(graph, edge.loc_to),
        edge.obs_mask,
    };
}

Mwpm& cluster_mwpm(ClusterSubgraph& cluster_subgraph, bool edge_correlations) {
    return cluster_subgraph_mwpm(cluster_subgraph, edge_correlations);
}


bool stopped_node_can_remap_to_target(
    const ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    size_t source_local_node_id) {
    if (source_local_node_id == NO_STOPPED_NODE) {
        return true;
    }
    if (target_subgraph.local_node_ids_are_global_ids && source_subgraph.local_node_ids_are_global_ids) {
        return source_local_node_id < target_subgraph.local_to_global_node_ids.size();
    }
    if (source_local_node_id >= source_subgraph.local_to_global_node_ids.size()) {
        return false;
    }
    size_t global_node_id = source_subgraph.local_to_global_node_ids[source_local_node_id];
    if (global_node_id >= target_subgraph.global_to_local_node_ids.size()) {
        return false;
    }
    size_t target_local_node_id = target_subgraph.global_to_local_node_ids[global_node_id];
    return target_local_node_id != SIZE_MAX;
}

bool stopped_edge_can_remap_to_target(
    const ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    const StoppedCompressedEdge& edge) {
    return stopped_node_can_remap_to_target(target_subgraph, source_subgraph, edge.loc_from) &&
        stopped_node_can_remap_to_target(target_subgraph, source_subgraph, edge.loc_to);
}

bool stopped_node_maps_to_unowned_target(
    ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    size_t source_local_node_id,
    bool edge_correlations) {
    if (source_local_node_id == NO_STOPPED_NODE) {
        return true;
    }
    size_t target_local_node_id = SIZE_MAX;
    if (target_subgraph.local_node_ids_are_global_ids && source_subgraph.local_node_ids_are_global_ids) {
        if (source_local_node_id >= target_subgraph.local_to_global_node_ids.size()) {
            return false;
        }
        target_local_node_id = source_local_node_id;
    } else {
        if (source_local_node_id >= source_subgraph.local_to_global_node_ids.size()) {
            return false;
        }
        size_t global_node_id = source_subgraph.local_to_global_node_ids[source_local_node_id];
        if (global_node_id >= target_subgraph.global_to_local_node_ids.size()) {
            return false;
        }
        target_local_node_id = target_subgraph.global_to_local_node_ids[global_node_id];
        if (target_local_node_id == SIZE_MAX) {
            return false;
        }
    }
    auto& target_mwpm = cluster_mwpm(target_subgraph, edge_correlations);
    if (target_local_node_id >= target_mwpm.flooder.graph.nodes.size()) {
        return false;
    }
    return target_mwpm.flooder.graph.nodes[target_local_node_id].region_that_arrived == nullptr;
}

size_t remap_stopped_node_to_target_local_id(
    const ClusterSubgraph& target_subgraph, const ClusterSubgraph& source_subgraph, size_t source_local_node_id) {
    if (source_local_node_id == NO_STOPPED_NODE) {
        return NO_STOPPED_NODE;
    }
    if (target_subgraph.local_node_ids_are_global_ids && source_subgraph.local_node_ids_are_global_ids) {
        if (source_local_node_id >= target_subgraph.local_to_global_node_ids.size()) {
            throw std::invalid_argument("Stopped-state node id is outside the full-graph target cluster subgraph.");
        }
        return source_local_node_id;
    }
    if (source_local_node_id >= source_subgraph.local_to_global_node_ids.size()) {
        throw std::invalid_argument("Stopped-state node id is outside the source cluster subgraph.");
    }
    size_t global_node_id = source_subgraph.local_to_global_node_ids[source_local_node_id];
    if (global_node_id >= target_subgraph.global_to_local_node_ids.size()) {
        throw std::invalid_argument("Stopped-state node id is outside the target cluster subgraph.");
    }
    size_t target_local_node_id = target_subgraph.global_to_local_node_ids[global_node_id];
    if (target_local_node_id == SIZE_MAX) {
        throw std::invalid_argument(
            "Stopped-state node is not present in the target cluster subgraph: target_cluster=" +
            std::to_string(target_subgraph.cluster_id) + " source_cluster=" + std::to_string(source_subgraph.cluster_id) +
            " source_local=" + std::to_string(source_local_node_id) + " global=" + std::to_string(global_node_id) +
            " target_local_nodes=" + std::to_string(target_subgraph.local_to_global_node_ids.size()) +
            " source_local_nodes=" + std::to_string(source_subgraph.local_to_global_node_ids.size()));
    }
    return target_local_node_id;
}

DetectorNode* remap_stopped_node_pointer(
    ClusterSubgraph& target_subgraph, const ClusterSubgraph& source_subgraph, size_t source_local_node_id,
    bool edge_correlations = false) {
    if (source_local_node_id == NO_STOPPED_NODE) {
        return nullptr;
    }
    size_t target_local_node_id = remap_stopped_node_to_target_local_id(
        target_subgraph, source_subgraph, source_local_node_id);
    return &cluster_mwpm(target_subgraph, edge_correlations).flooder.graph.nodes[target_local_node_id];
}

DetectorNode* remap_stopped_node_pointer(
    const ClusterSubgraph& target_subgraph, const ClusterSubgraph& source_subgraph, size_t source_local_node_id) {
    return remap_stopped_node_pointer(
        const_cast<ClusterSubgraph&>(target_subgraph), source_subgraph, source_local_node_id, false);
}

CompressedEdge import_stopped_compressed_edge(
    const ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    const StoppedCompressedEdge& edge) {
    return {
        remap_stopped_node_pointer(target_subgraph, source_subgraph, edge.loc_from),
        remap_stopped_node_pointer(target_subgraph, source_subgraph, edge.loc_to),
        edge.obs_mask,
    };
}


MwpmStoppedState filter_stopped_state_to_target_buffer(
    ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    const MwpmStoppedState& source_state,
    bool edge_correlations) {
    MwpmStoppedState filtered;
    filtered.algorithmic_time = source_state.algorithmic_time;
    filtered.num_nodes = source_state.num_nodes;
    filtered.num_observables = source_state.num_observables;
    filtered.negative_weight_detection_events = source_state.negative_weight_detection_events;
    filtered.negative_weight_observables = source_state.negative_weight_observables;
    filtered.negative_weight_obs_mask = source_state.negative_weight_obs_mask;
    filtered.negative_weight_sum = source_state.negative_weight_sum;

    std::vector<uint8_t> keep_region(source_state.region_states.size(), 1);
    for (const auto& region_state : source_state.region_states) {
        if (region_state.region_id >= keep_region.size()) {
            continue;
        }
        bool keep = true;
        if (!stopped_edge_can_remap_to_target(target_subgraph, source_subgraph, region_state.match.edge)) {
            keep = false;
        }
        for (auto source_local_node_id : region_state.shell_area_node_ids) {
            if (!stopped_node_can_remap_to_target(target_subgraph, source_subgraph, source_local_node_id)) {
                keep = false;
                break;
            }
        }
        if (keep) {
            for (const auto& child : region_state.blossom_children) {
                if (!stopped_edge_can_remap_to_target(target_subgraph, source_subgraph, child.edge)) {
                    keep = false;
                    break;
                }
            }
        }
        keep_region[region_state.region_id] = keep ? 1 : 0;
    }

    // If importing a detector state would overwrite existing target ownership,
    // drop the whole associated region.  This makes repeated broadcasts and
    // descendant/ancestor duplicate paths naturally become no-ops instead of
    // triggering the parent-overwrite recovery path.
    for (const auto& detector_state : source_state.detector_states) {
        if (detector_state.node_id == NO_STOPPED_NODE || detector_state.region_id == NO_STOPPED_REGION) {
            continue;
        }
        if (detector_state.region_id >= keep_region.size()) {
            continue;
        }
        if (!stopped_node_maps_to_unowned_target(
                target_subgraph, source_subgraph, detector_state.node_id, edge_correlations)) {
            keep_region[detector_state.region_id] = 0;
        }
    }

    bool changed = true;
    while (changed) {
        changed = false;
        for (const auto& region_state : source_state.region_states) {
            if (region_state.region_id >= keep_region.size() || !keep_region[region_state.region_id]) {
                continue;
            }
            auto require_region = [&](size_t other_id) {
                if (other_id == NO_STOPPED_REGION) {
                    return true;
                }
                return other_id < keep_region.size() && keep_region[other_id];
            };
            bool keep = true;
            keep &= require_region(region_state.blossom_parent_id);
            keep &= require_region(region_state.blossom_parent_top_id);
            keep &= require_region(region_state.match.region_id);
            for (const auto& child : region_state.blossom_children) {
                keep &= require_region(child.region_id);
            }
            if (!keep) {
                keep_region[region_state.region_id] = 0;
                changed = true;
            }
        }
    }

    std::vector<size_t> region_id_map(source_state.region_states.size(), NO_STOPPED_REGION);
    for (const auto& region_state : source_state.region_states) {
        if (region_state.region_id >= keep_region.size() || !keep_region[region_state.region_id]) {
            continue;
        }
        region_id_map[region_state.region_id] = filtered.region_states.size();
        filtered.region_states.push_back(region_state);
    }
    filtered.num_regions = filtered.region_states.size();

    auto map_region_id = [&](size_t old_id) {
        if (old_id == NO_STOPPED_REGION) {
            return NO_STOPPED_REGION;
        }
        if (old_id >= region_id_map.size()) {
            return NO_STOPPED_REGION;
        }
        return region_id_map[old_id];
    };
    for (auto& region_state : filtered.region_states) {
        region_state.region_id = map_region_id(region_state.region_id);
        region_state.blossom_parent_id = map_region_id(region_state.blossom_parent_id);
        region_state.blossom_parent_top_id = map_region_id(region_state.blossom_parent_top_id);
        region_state.match.region_id = map_region_id(region_state.match.region_id);
        for (auto& child : region_state.blossom_children) {
            child.region_id = map_region_id(child.region_id);
        }
    }

    for (auto detector_state : source_state.detector_states) {
        if (detector_state.node_id == NO_STOPPED_NODE || detector_state.region_id == NO_STOPPED_REGION) {
            continue;
        }
        if (!stopped_node_maps_to_unowned_target(
                target_subgraph, source_subgraph, detector_state.node_id, edge_correlations)) {
            continue;
        }
        detector_state.region_id = map_region_id(detector_state.region_id);
        detector_state.top_region_id = map_region_id(detector_state.top_region_id);
        if (detector_state.region_id == NO_STOPPED_REGION) {
            continue;
        }
        filtered.detector_states.push_back(detector_state);
    }

    return filtered;
}


MwpmStoppedState filter_stopped_state_to_reached_sources(
    const MwpmStoppedState& source_state,
    const std::vector<uint64_t>& active_detector_local_node_ids) {
    MwpmStoppedState filtered;
    filtered.algorithmic_time = source_state.algorithmic_time;
    filtered.num_nodes = source_state.num_nodes;
    filtered.num_observables = source_state.num_observables;
    filtered.negative_weight_detection_events = source_state.negative_weight_detection_events;
    filtered.negative_weight_observables = source_state.negative_weight_observables;
    filtered.negative_weight_obs_mask = source_state.negative_weight_obs_mask;
    filtered.negative_weight_sum = source_state.negative_weight_sum;

    if (active_detector_local_node_ids.empty() || source_state.region_states.empty()) {
        return filtered;
    }

    std::unordered_set<size_t> active_sources;
    active_sources.reserve(active_detector_local_node_ids.size() * 2 + 1);
    for (auto node_id : active_detector_local_node_ids) {
        active_sources.insert(static_cast<size_t>(node_id));
    }

    std::vector<uint8_t> keep_detector(source_state.detector_states.size(), 0);
    std::vector<uint8_t> keep_region(source_state.region_states.size(), 0);
    std::unordered_set<size_t> kept_detector_nodes;
    kept_detector_nodes.reserve(source_state.detector_states.size());

    auto mark_region = [&](size_t region_id) {
        if (region_id != NO_STOPPED_REGION && region_id < keep_region.size()) {
            keep_region[region_id] = 1;
        }
    };

    for (size_t k = 0; k < source_state.detector_states.size(); k++) {
        const auto& detector_state = source_state.detector_states[k];
        if (detector_state.node_id == NO_STOPPED_NODE || detector_state.region_id == NO_STOPPED_REGION) {
            continue;
        }
        if (!active_sources.contains(detector_state.reached_from_source_id)) {
            continue;
        }
        keep_detector[k] = 1;
        kept_detector_nodes.insert(detector_state.node_id);
        mark_region(detector_state.region_id);
        mark_region(detector_state.top_region_id);
    }

    // Keep the structural closure of each retained region.  Blossom children and
    // ancestry are local representation details that must stay internally
    // consistent after compacting region ids.  Matches are not used to pull in a
    // new active-source component: if the partner region was not itself reached
    // from this cluster's active detectors, the match is cleared below instead of
    // exporting unrelated state.
    bool changed = true;
    while (changed) {
        changed = false;
        for (const auto& region_state : source_state.region_states) {
            if (region_state.region_id >= keep_region.size() || !keep_region[region_state.region_id]) {
                continue;
            }
            auto add_dependency = [&](size_t other_id) {
                if (other_id == NO_STOPPED_REGION || other_id >= keep_region.size() || keep_region[other_id]) {
                    return;
                }
                keep_region[other_id] = 1;
                changed = true;
            };
            add_dependency(region_state.blossom_parent_id);
            add_dependency(region_state.blossom_parent_top_id);
            for (const auto& child : region_state.blossom_children) {
                add_dependency(child.region_id);
            }
        }
    }

    std::vector<size_t> region_id_map(source_state.region_states.size(), NO_STOPPED_REGION);
    for (const auto& region_state : source_state.region_states) {
        if (region_state.region_id >= keep_region.size() || !keep_region[region_state.region_id]) {
            continue;
        }
        region_id_map[region_state.region_id] = filtered.region_states.size();
        filtered.region_states.push_back(region_state);
    }
    filtered.num_regions = filtered.region_states.size();

    auto map_region_id = [&](size_t old_id) {
        if (old_id == NO_STOPPED_REGION) {
            return NO_STOPPED_REGION;
        }
        if (old_id >= region_id_map.size()) {
            return NO_STOPPED_REGION;
        }
        return region_id_map[old_id];
    };

    for (auto& region_state : filtered.region_states) {
        region_state.region_id = map_region_id(region_state.region_id);
        region_state.blossom_parent_id = map_region_id(region_state.blossom_parent_id);
        region_state.blossom_parent_top_id = map_region_id(region_state.blossom_parent_top_id);

        auto mapped_match_region_id = map_region_id(region_state.match.region_id);
        if (region_state.match.region_id != NO_STOPPED_REGION && mapped_match_region_id == NO_STOPPED_REGION) {
            region_state.match = StoppedMatch{};
        } else {
            region_state.match.region_id = mapped_match_region_id;
        }

        std::vector<StoppedRegionEdge> kept_children;
        kept_children.reserve(region_state.blossom_children.size());
        for (auto child : region_state.blossom_children) {
            child.region_id = map_region_id(child.region_id);
            if (child.region_id != NO_STOPPED_REGION) {
                kept_children.push_back(child);
            }
        }
        region_state.blossom_children = std::move(kept_children);

        std::vector<size_t> kept_shell;
        kept_shell.reserve(region_state.shell_area_node_ids.size());
        for (auto node_id : region_state.shell_area_node_ids) {
            if (kept_detector_nodes.contains(node_id)) {
                kept_shell.push_back(node_id);
            }
        }
        region_state.shell_area_node_ids = std::move(kept_shell);
    }

    for (size_t k = 0; k < source_state.detector_states.size(); k++) {
        if (!keep_detector[k]) {
            continue;
        }
        auto detector_state = source_state.detector_states[k];
        detector_state.region_id = map_region_id(detector_state.region_id);
        detector_state.top_region_id = map_region_id(detector_state.top_region_id);
        if (detector_state.region_id == NO_STOPPED_REGION) {
            continue;
        }
        filtered.detector_states.push_back(detector_state);
    }

    if (filtered.detector_states.empty()) {
        filtered.region_states.clear();
        filtered.num_regions = 0;
    }
    return filtered;
}

void validate_stopped_state_import_does_not_overwrite_detector_ownership(
    const ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    const MwpmStoppedState& source_state) {
    for (const auto& detector_state : source_state.detector_states) {
        if (detector_state.node_id == NO_STOPPED_NODE || detector_state.region_id == NO_STOPPED_REGION) {
            continue;
        }
        auto* target_node = remap_stopped_node_pointer(target_subgraph, source_subgraph, detector_state.node_id);
        if (target_node->region_that_arrived != nullptr) {
            throw std::invalid_argument("Stopped-state import would overwrite detector ownership in the target cluster.");
        }
    }
}

void validate_quiescent_state(const Mwpm& mwpm) {
    // The flooder queue can contain stale invalid tentative events even after
    // sparse blossom is quiescent.  The stopped state only requires that no
    // valid tentative events and no live alternating-tree nodes remain.
    if (mwpm.flooder.has_valid_tentative_events()) {
        throw std::invalid_argument("Cannot export Mwpm stopped state while valid flooder events remain.");
    }

    auto live_alt_tree_nodes = collect_live_arena_objects(mwpm.node_arena);
    if (!live_alt_tree_nodes.empty()) {
        throw std::invalid_argument(
            "Cannot export Mwpm stopped state while alternating-tree structure is still live.");
    }
}

}  // namespace

MwpmStoppedState export_mwpm_stopped_state(const Mwpm& mwpm) {
    validate_quiescent_state(mwpm);

    MwpmStoppedState state;
    state.algorithmic_time = mwpm.flooder.queue.cur_time;
    state.num_nodes = mwpm.flooder.graph.nodes.size();
    state.num_observables = mwpm.flooder.graph.num_observables;
    state.negative_weight_detection_events = mwpm.flooder.negative_weight_detection_events;
    state.negative_weight_observables = mwpm.flooder.negative_weight_observables;
    state.negative_weight_obs_mask = mwpm.flooder.negative_weight_obs_mask;
    state.negative_weight_sum = mwpm.flooder.negative_weight_sum;

    auto live_regions = collect_live_arena_objects(mwpm.flooder.region_arena);
    state.num_regions = live_regions.size();

    std::unordered_map<const GraphFillRegion*, size_t> region_ids;
    region_ids.reserve(live_regions.size());
    for (size_t k = 0; k < live_regions.size(); k++) {
        region_ids.emplace(live_regions[k], k);
    }

    // Sparse stopped-state export: only nodes actually owned by a stopped region
    // are needed for child-to-parent import. Empty detector nodes have no flooder
    // ownership to restore, and exporting them made each child import scale with
    // the full worker graph size in full-graph-worker mode.
    state.detector_states.reserve(live_regions.size() * 4);
    for (size_t node_id = 0; node_id < mwpm.flooder.graph.nodes.size(); node_id++) {
        const auto& node = mwpm.flooder.graph.nodes[node_id];
        if (node.region_that_arrived == nullptr) {
            continue;
        }
        try {
            state.detector_states.push_back({
                node_id,
                stopped_region_id(node.region_that_arrived, region_ids),
                stopped_region_id(node.region_that_arrived_top, region_ids),
                stopped_node_id(mwpm.flooder.graph, node.reached_from_source),
                node.observables_crossed_from_source,
                node.radius_of_arrival,
                node.wrapped_radius_cached,
            });
        } catch (const std::invalid_argument& ex) {
            throw std::invalid_argument(
                "Failed exporting DetectorNodeStoppedState for node " + std::to_string(node_id) + ": " + ex.what());
        }
    }

    state.region_states.reserve(live_regions.size());
    for (size_t region_id = 0; region_id < live_regions.size(); region_id++) {
        const auto* region = live_regions[region_id];
        if (region->alt_tree_node != nullptr) {
            throw std::invalid_argument(
                "Cannot export Mwpm stopped state while a region still belongs to an alternating tree.");
        }
        GraphFillRegionStoppedState region_state;
        region_state.region_id = region_id;
        try {
            region_state.blossom_parent_id = stopped_region_id(region->blossom_parent, region_ids);
            region_state.blossom_parent_top_id = stopped_region_id(region->blossom_parent_top, region_ids);
        } catch (const std::invalid_argument& ex) {
            throw std::invalid_argument(
                "Failed exporting blossom ancestry for region " + std::to_string(region_id) + ": " + ex.what());
        }
        region_state.alt_tree_node_id = NO_STOPPED_ALT_TREE_NODE;
        region_state.radius = region->radius;
        if (region->blossom_parent == nullptr) {
            try {
                region_state.match = {
                    stopped_region_id(region->match.region, region_ids),
                    export_stopped_compressed_edge(mwpm.flooder.graph, region->match.edge),
                };
            } catch (const std::invalid_argument& ex) {
                throw std::invalid_argument(
                    "Failed exporting match for top-level region " + std::to_string(region_id) + ": " + ex.what());
            }
        }
        region_state.blossom_children.reserve(region->blossom_children.size());
        for (const auto& child : region->blossom_children) {
            try {
                region_state.blossom_children.push_back({
                    stopped_region_id(child.region, region_ids),
                    export_stopped_compressed_edge(mwpm.flooder.graph, child.edge),
                });
            } catch (const std::invalid_argument& ex) {
                throw std::invalid_argument(
                    "Failed exporting blossom child for region " + std::to_string(region_id) + ": " + ex.what());
            }
        }
        region_state.shell_area_node_ids.reserve(region->shell_area.size());
        for (const auto* node : region->shell_area) {
            region_state.shell_area_node_ids.push_back(stopped_node_id(mwpm.flooder.graph, node));
        }
        state.region_states.push_back(std::move(region_state));
    }

    return state;
}


MwpmStoppedState export_mwpm_stopped_state_for_active_detectors(
    const Mwpm& mwpm,
    const std::vector<uint64_t>& active_detector_local_node_ids) {
    auto state = export_mwpm_stopped_state(mwpm);
    auto filtered = filter_stopped_state_to_reached_sources(state, active_detector_local_node_ids);
    if (std::getenv("PYMATCHING_PARENTLESS_EXPORT_STATS") != nullptr) {
        std::cerr << "PARENTLESS_ACTIVE_EXPORT active_sources=" << active_detector_local_node_ids.size()
                  << " detector_states_before=" << state.detector_states.size()
                  << " detector_states_after=" << filtered.detector_states.size()
                  << " region_states_before=" << state.region_states.size()
                  << " region_states_after=" << filtered.region_states.size() << "\n";
    }
    return filtered;
}

MwpmStoppedState decode_detection_events_to_stopped_state(
    Mwpm& mwpm, const std::vector<uint64_t>& detection_events, bool edge_correlations) {
    if (edge_correlations) {
        throw std::invalid_argument(
            "Stopped-state decoding does not yet support edge correlations.");
    }
    pm::process_timeline_until_completion(mwpm, detection_events);
    return export_mwpm_stopped_state(mwpm);
}

std::vector<uint64_t> remap_source_detection_events_from_stopped_state(
    const ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    const MwpmStoppedState& source_state) {
    std::vector<uint64_t> remapped;
    remapped.reserve(source_state.detector_states.size());
    for (const auto& detector_state : source_state.detector_states) {
        if (detector_state.region_id == NO_STOPPED_REGION) {
            continue;
        }
        if (detector_state.reached_from_source_id != detector_state.node_id) {
            continue;
        }
        auto* target_node = remap_stopped_node_pointer(target_subgraph, source_subgraph, detector_state.node_id);
        auto target_node_id = static_cast<uint64_t>(
            target_node - cluster_mwpm(const_cast<ClusterSubgraph&>(target_subgraph), false).flooder.graph.nodes.data());
        remapped.push_back(target_node_id);
    }
    std::sort(remapped.begin(), remapped.end());
    remapped.erase(std::unique(remapped.begin(), remapped.end()), remapped.end());
    return remapped;
}

bool stopped_state_import_would_overwrite_detector_ownership(
    const ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    const MwpmStoppedState& source_state) {
    try {
        validate_stopped_state_import_does_not_overwrite_detector_ownership(
            target_subgraph, source_subgraph, source_state);
        return false;
    } catch (const std::invalid_argument&) {
        return true;
    }
}

void import_stopped_state_into_cluster_subgraph(
    ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    const MwpmStoppedState& source_state,
    bool edge_correlations) {
    import_stopped_state_into_cluster_subgraph(
        target_subgraph,
        source_subgraph,
        source_state,
        StoppedStateImportOptions{},
        edge_correlations);
}

void import_stopped_state_into_cluster_subgraph(
    ClusterSubgraph& target_subgraph,
    const ClusterSubgraph& source_subgraph,
    const MwpmStoppedState& source_state,
    const StoppedStateImportOptions& options,
    bool edge_correlations) {
    if (edge_correlations) {
        throw std::invalid_argument("Stopped-state import does not yet support edge correlations.");
    }

    MwpmStoppedState filtered_state_storage;
    const MwpmStoppedState* import_state = &source_state;
    if (options.skip_unmappable_or_conflicting_state) {
        filtered_state_storage = filter_stopped_state_to_target_buffer(
            target_subgraph, source_subgraph, source_state, edge_correlations);
        import_state = &filtered_state_storage;
        if (import_state->detector_states.empty() || import_state->region_states.empty()) {
            return;
        }
    }

    auto& target_mwpm = cluster_mwpm(target_subgraph, edge_correlations);
    // Normal parent-child imports reject detector-ownership overwrites because they
    // indicate the parent has already evolved through the child's support.  The
    // experimental parentless broadcast path intentionally permits duplicate
    // active-detector-limited states to arrive through multiple upper levels; in
    // that mode the later import overwrites the detector ownership.
    if (!options.allow_detector_ownership_overwrite) {
        validate_stopped_state_import_does_not_overwrite_detector_ownership(
            target_subgraph, source_subgraph, *import_state);
    }

    // Only advance the parent flooder clock when the parent queue is empty.
    // The radix heap buckets are keyed relative to queue.cur_time; mutating
    // cur_time while events are still queued invalidates the bucket layout and
    // can later make dequeue walk past the bucket array.  If the parent has
    // queued work, the lockstep scheduler will advance the queue through its
    // normal dequeue path.
    if (target_mwpm.flooder.queue.empty()) {
        target_mwpm.flooder.queue.cur_time = std::max(
            target_mwpm.flooder.queue.cur_time, import_state->algorithmic_time);
    }

    std::vector<GraphFillRegion*> imported_regions(import_state->region_states.size(), nullptr);
    for (size_t region_id = 0; region_id < import_state->region_states.size(); region_id++) {
        auto* region = target_mwpm.flooder.region_arena.alloc_default_constructed();
        region->alt_tree_node = nullptr;
        region->radius = import_state->region_states[region_id].radius;
        region->match.clear();
        region->blossom_children.clear();
        region->shell_area.clear();
        imported_regions[region_id] = region;
    }

    for (const auto& region_state : import_state->region_states) {
        auto* region = imported_regions[region_state.region_id];
        region->blossom_parent = region_state.blossom_parent_id == NO_STOPPED_REGION
                                     ? nullptr
                                     : imported_regions[region_state.blossom_parent_id];
        region->blossom_parent_top = region_state.blossom_parent_top_id == NO_STOPPED_REGION
                                         ? region
                                         : imported_regions[region_state.blossom_parent_top_id];
    }

    for (const auto& region_state : import_state->region_states) {
        auto* region = imported_regions[region_state.region_id];
        if (region_state.match.region_id != NO_STOPPED_REGION || region_state.match.edge.loc_from != NO_STOPPED_NODE ||
            region_state.match.edge.loc_to != NO_STOPPED_NODE) {
            region->match = {
                region_state.match.region_id == NO_STOPPED_REGION ? nullptr : imported_regions[region_state.match.region_id],
                import_stopped_compressed_edge(target_subgraph, source_subgraph, region_state.match.edge),
            };
        }
        for (const auto& child_state : region_state.blossom_children) {
            region->blossom_children.push_back({
                imported_regions[child_state.region_id],
                import_stopped_compressed_edge(target_subgraph, source_subgraph, child_state.edge),
            });
        }
        for (auto source_local_node_id : region_state.shell_area_node_ids) {
            auto* node = remap_stopped_node_pointer(target_subgraph, source_subgraph, source_local_node_id, edge_correlations);
            region->shell_area.push_back(node);
        }
    }

    for (const auto& detector_state : import_state->detector_states) {
        if (detector_state.node_id == NO_STOPPED_NODE || detector_state.region_id == NO_STOPPED_REGION) {
            continue;
        }
        auto target_local_node_id = remap_stopped_node_to_target_local_id(
            target_subgraph, source_subgraph, detector_state.node_id);
        auto* target_node = &target_mwpm.flooder.graph.nodes[target_local_node_id];
        target_node->region_that_arrived = imported_regions[detector_state.region_id];
        target_node->region_that_arrived_top = detector_state.top_region_id == NO_STOPPED_REGION
                                                   ? nullptr
                                                   : imported_regions[detector_state.top_region_id];
        target_node->reached_from_source = remap_stopped_node_pointer(
            target_subgraph, source_subgraph, detector_state.reached_from_source_id, edge_correlations);
        target_node->observables_crossed_from_source = detector_state.observables_crossed_from_source;
        target_node->radius_of_arrival = detector_state.radius_of_arrival;
        target_node->wrapped_radius_cached = detector_state.wrapped_radius_cached;
        target_node->node_event_tracker.clear();
        if (options.imported_target_local_node_ids != nullptr) {
            options.imported_target_local_node_ids->push_back(target_local_node_id);
        }
    }

    if (options.reschedule_after_import) {
        std::vector<size_t> target_local_node_ids;
        target_local_node_ids.reserve(import_state->detector_states.size());
        if (options.imported_target_local_node_ids != nullptr) {
            target_local_node_ids = *options.imported_target_local_node_ids;
        } else {
            for (const auto& detector_state : import_state->detector_states) {
                if (detector_state.node_id == NO_STOPPED_NODE || detector_state.region_id == NO_STOPPED_REGION) {
                    continue;
                }
                target_local_node_ids.push_back(remap_stopped_node_to_target_local_id(
                    target_subgraph, source_subgraph, detector_state.node_id));
            }
        }
        reschedule_imported_detector_nodes(target_subgraph, target_local_node_ids, edge_correlations);
    }
}

void reschedule_imported_detector_nodes(
    ClusterSubgraph& target_subgraph,
    const std::vector<size_t>& target_local_node_ids,
    bool edge_correlations) {
    if (edge_correlations) {
        throw std::invalid_argument("Stopped-state import reschedule does not yet support edge correlations.");
    }
    if (target_local_node_ids.empty()) {
        return;
    }
    auto& target_mwpm = cluster_mwpm(target_subgraph, edge_correlations);
    std::vector<size_t> unique_ids = target_local_node_ids;
    std::sort(unique_ids.begin(), unique_ids.end());
    unique_ids.erase(std::unique(unique_ids.begin(), unique_ids.end()), unique_ids.end());
    for (auto target_local_node_id : unique_ids) {
        if (target_local_node_id >= target_mwpm.flooder.graph.nodes.size()) {
            throw std::invalid_argument("Imported detector id is outside target MWPM graph.");
        }
        target_mwpm.flooder.reschedule_events_at_detector_node(
            target_mwpm.flooder.graph.nodes[target_local_node_id]);
    }
}

}  // namespace pm
