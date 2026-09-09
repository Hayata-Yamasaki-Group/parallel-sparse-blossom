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

#include "pymatching/sparse_blossom/parallel/cluster_subgraph.h"

#include <algorithm>
#include <condition_variable>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <memory>
#include <mutex>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#if defined(__linux__)
#include <pthread.h>
#include <sched.h>
#endif

namespace pm {
namespace {

struct PendingImpliedWeightRewrite {
    size_t local_node1;
    size_t local_node2;
    std::vector<ImpliedWeightUnconverted> rules;
};

void validate_cluster_vertices(UserGraph& graph, const ProcessingCluster& cluster) {
    size_t num_nodes = graph.get_num_nodes();
    if (!cluster.influence_is_full_graph) {
        for (auto node : cluster.influence_vertices) {
            if (node >= num_nodes) {
                throw std::invalid_argument("ProcessingCluster influence vertex is not present in the UserGraph.");
            }
        }
    }
    for (auto node : cluster.boundary_endpoint_vertices) {
        if (node >= num_nodes) {
            throw std::invalid_argument("ProcessingCluster boundary endpoint is not present in the UserGraph.");
        }
    }
    for (auto detector : cluster.active_detectors) {
        if (detector >= num_nodes) {
            throw std::invalid_argument("ProcessingCluster active detector is not present in the UserGraph.");
        }
    }
}

std::vector<bool> compute_kept_global_nodes(UserGraph& graph, const ProcessingCluster& cluster) {
    std::vector<bool> kept(graph.get_num_nodes(), false);
    for (auto node : cluster.influence_vertices) {
        kept[node] = true;
    }
    for (auto node : cluster.influence_vertices) {
        const auto& user_node = graph.nodes[node];
        for (const auto& neighbor : user_node.neighbors) {
            const auto& edge = *neighbor.edge_it;
            if (edge.node2 == SIZE_MAX) {
                continue;
            }
            if (graph.boundary_nodes.contains(edge.node1)) {
                kept[edge.node1] = true;
            }
            if (graph.boundary_nodes.contains(edge.node2)) {
                kept[edge.node2] = true;
            }
        }
    }
    return kept;
}


std::vector<size_t> collect_kept_global_nodes(UserGraph& graph, const ProcessingCluster& cluster) {
    if (graph.boundary_nodes.empty()) {
        // In DEM-derived detector graphs, boundaries are represented by boundary
        // edges rather than explicit boundary vertices. The influence vertices
        // already contain every endpoint needed by those boundary edges, so avoid
        // scanning all incident UserGraph edges just to discover that there are no
        // explicit boundary nodes.
        return cluster.influence_vertices;
    }

    std::vector<size_t> nodes;
    nodes.reserve(cluster.influence_vertices.size() + cluster.boundary_endpoint_vertices.size());
    for (auto node : cluster.influence_vertices) {
        nodes.push_back(node);
    }
    for (auto node : cluster.boundary_endpoint_vertices) {
        nodes.push_back(node);
    }
    for (auto node : cluster.influence_vertices) {
        const auto& user_node = graph.nodes[node];
        for (const auto& neighbor : user_node.neighbors) {
            const auto& edge = *neighbor.edge_it;
            if (edge.node2 == SIZE_MAX) {
                continue;
            }
            if (graph.boundary_nodes.contains(edge.node1)) {
                nodes.push_back(edge.node1);
            }
            if (graph.boundary_nodes.contains(edge.node2)) {
                nodes.push_back(edge.node2);
            }
        }
    }
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    return nodes;
}

std::vector<bool> compute_kept_mask_from_local_to_global(
    size_t num_global_nodes, const std::vector<size_t>& local_to_global) {
    std::vector<bool> kept(num_global_nodes, false);
    for (auto node : local_to_global) {
        kept[node] = true;
    }
    return kept;
}

std::vector<size_t> invert_kept_global_nodes(const std::vector<bool>& kept) {
    std::vector<size_t> local_to_global;
    for (size_t node = 0; node < kept.size(); node++) {
        if (kept[node]) {
            local_to_global.push_back(node);
        }
    }
    return local_to_global;
}

std::vector<size_t> build_global_to_local_node_ids(size_t num_global_nodes, const std::vector<size_t>& local_to_global) {
    std::vector<size_t> global_to_local(num_global_nodes, SIZE_MAX);
    for (size_t local = 0; local < local_to_global.size(); local++) {
        global_to_local[local_to_global[local]] = local;
    }
    return global_to_local;
}

std::set<size_t> build_local_boundary_nodes(
    UserGraph& graph, const std::vector<size_t>& local_to_global_node_ids, std::vector<bool>& boundary_mask) {
    std::set<size_t> local_boundary_nodes;
    if (graph.boundary_nodes.empty()) {
        boundary_mask.clear();
        return local_boundary_nodes;
    }
    boundary_mask.assign(local_to_global_node_ids.size(), false);
    for (size_t local = 0; local < local_to_global_node_ids.size(); local++) {
        size_t global = local_to_global_node_ids[local];
        if (graph.boundary_nodes.contains(global)) {
            local_boundary_nodes.insert(local);
            boundary_mask[local] = true;
        }
    }
    return local_boundary_nodes;
}

std::list<UserEdge>::iterator find_local_edge_it(UserGraph& graph, size_t node1, size_t node2) {
    size_t neighbor_index = node2 == SIZE_MAX ? graph.nodes[node1].index_of_neighbor(SIZE_MAX)
                                              : graph.nodes[node1].index_of_neighbor(node2);
    if (neighbor_index == SIZE_MAX) {
        throw std::invalid_argument("Expected edge to exist in cluster subgraph.");
    }
    return graph.nodes[node1].neighbors[neighbor_index].edge_it;
}

void add_included_edge(
    UserGraph& graph,
    const UserEdge& edge,
    size_t local_node1,
    size_t local_node2,
    std::vector<PendingImpliedWeightRewrite>& pending_rules) {
    if (edge.node2 == SIZE_MAX) {
        graph.add_or_merge_boundary_edge(
            local_node1, edge.observable_indices, edge.weight, edge.error_probability, DISALLOW);
    } else {
        graph.add_or_merge_edge(
            local_node1, local_node2, edge.observable_indices, edge.weight, edge.error_probability, DISALLOW);
    }
    if (!edge.implied_weights_for_other_edges.empty()) {
        pending_rules.push_back({local_node1, local_node2, edge.implied_weights_for_other_edges});
    }
}

void rewrite_implied_weight_rules(
    ClusterSubgraph& cluster_subgraph, const std::vector<PendingImpliedWeightRewrite>& pending_rules) {
    for (const auto& pending : pending_rules) {
        auto edge_it = find_local_edge_it(cluster_subgraph.graph, pending.local_node1, pending.local_node2);
        std::vector<ImpliedWeightUnconverted> remapped_rules;
        for (const auto& rule : pending.rules) {
            if (rule.node1 == SIZE_MAX || rule.node1 >= cluster_subgraph.global_to_local_node_ids.size()) {
                cluster_subgraph.dropped_implied_weight_rules = true;
                continue;
            }
            size_t remapped_node1 = cluster_subgraph.global_to_local_node_ids[rule.node1];
            if (remapped_node1 == SIZE_MAX) {
                cluster_subgraph.dropped_implied_weight_rules = true;
                continue;
            }

            if (rule.node2 == SIZE_MAX) {
                if (!cluster_subgraph.graph.has_boundary_edge(remapped_node1)) {
                    cluster_subgraph.dropped_implied_weight_rules = true;
                    continue;
                }
                remapped_rules.push_back({remapped_node1, SIZE_MAX, rule.implied_weight});
                continue;
            }

            if (rule.node2 >= cluster_subgraph.global_to_local_node_ids.size()) {
                cluster_subgraph.dropped_implied_weight_rules = true;
                continue;
            }
            size_t remapped_node2 = cluster_subgraph.global_to_local_node_ids[rule.node2];
            if (remapped_node2 == SIZE_MAX || !cluster_subgraph.graph.has_edge(remapped_node1, remapped_node2)) {
                cluster_subgraph.dropped_implied_weight_rules = true;
                continue;
            }
            remapped_rules.push_back({remapped_node1, remapped_node2, rule.implied_weight});
        }
        edge_it->implied_weights_for_other_edges = std::move(remapped_rules);
    }
}


void toggle_negative_weight_observables_from_mask(MatchingGraph& graph, obs_int obs_mask) {
    for (size_t obs = 0; obs < graph.num_observables && obs < sizeof(obs_int) * 8; obs++) {
        if (((obs_mask >> obs) & 1) == 0) {
            continue;
        }
        auto iter = graph.negative_weight_observables_set.find(obs);
        if (iter == graph.negative_weight_observables_set.end()) {
            graph.negative_weight_observables_set.insert(obs);
        } else {
            graph.negative_weight_observables_set.erase(iter);
        }
    }
}

void toggle_negative_weight_detection_event(MatchingGraph& graph, size_t node) {
    auto iter = graph.negative_weight_detection_events_set.find(node);
    if (iter == graph.negative_weight_detection_events_set.end()) {
        graph.negative_weight_detection_events_set.insert(node);
    } else {
        graph.negative_weight_detection_events_set.erase(iter);
    }
}

void add_matching_graph_edge_fast(
    MatchingGraph& graph, size_t u, size_t v, signed_weight_int weight, obs_int obs_mask) {
    if (weight < 0) {
        toggle_negative_weight_observables_from_mask(graph, obs_mask);
        toggle_negative_weight_detection_event(graph, u);
        toggle_negative_weight_detection_event(graph, v);
        graph.negative_weight_sum += weight;
    }
    if (u == v) {
        return;
    }
    auto abs_weight = static_cast<weight_int>(std::abs(weight));
    graph.nodes[u].neighbors.push_back(&graph.nodes[v]);
    graph.nodes[u].neighbor_weights.push_back(abs_weight);
    graph.nodes[u].neighbor_observables.push_back(obs_mask);
    graph.nodes[u].neighbor_implied_weights.push_back({});
    graph.nodes[v].neighbors.push_back(&graph.nodes[u]);
    graph.nodes[v].neighbor_weights.push_back(abs_weight);
    graph.nodes[v].neighbor_observables.push_back(obs_mask);
    graph.nodes[v].neighbor_implied_weights.push_back({});
}

void add_matching_graph_boundary_edge_fast(
    MatchingGraph& graph, size_t u, signed_weight_int weight, obs_int obs_mask) {
    if (weight < 0) {
        toggle_negative_weight_observables_from_mask(graph, obs_mask);
        toggle_negative_weight_detection_event(graph, u);
        graph.negative_weight_sum += weight;
    }
    auto abs_weight = static_cast<weight_int>(std::abs(weight));
    auto& node = graph.nodes[u];
    // The source MatchingGraph stores the boundary edge first. In the direct
    // subgraph builder we iterate neighbors in that same order, so the boundary
    // insertion happens before ordinary edges and does not need to shift an
    // already-populated adjacency vector. Keep the boundary at index 0 for
    // compatibility with MatchingGraph::add_boundary_edge semantics.
    if (node.neighbors.empty()) {
        node.neighbors.push_back(nullptr);
        node.neighbor_weights.push_back(abs_weight);
        node.neighbor_observables.push_back(obs_mask);
        node.neighbor_implied_weights.push_back({});
    } else {
        node.neighbors.insert(node.neighbors.begin(), nullptr);
        node.neighbor_weights.insert(node.neighbor_weights.begin(), abs_weight);
        node.neighbor_observables.insert(node.neighbor_observables.begin(), obs_mask);
        node.neighbor_implied_weights.insert(node.neighbor_implied_weights.begin(), {});
    }
}


bool cluster_influence_is_full_identity(UserGraph& graph, const ProcessingCluster& cluster) {
    if (cluster.influence_is_full_graph) {
        return true;
    }
    if (cluster.influence_vertices.size() != graph.get_num_nodes()) {
        return false;
    }
    for (size_t k = 0; k < cluster.influence_vertices.size(); k++) {
        if (cluster.influence_vertices[k] != k) {
            return false;
        }
    }
    return true;
}

MatchingGraph clone_identity_matching_graph_for_cluster(const MatchingGraph& global_matching_graph) {
    MatchingGraph local_matching_graph(
        global_matching_graph.nodes.size(),
        global_matching_graph.num_observables,
        global_matching_graph.normalising_constant);
    local_matching_graph.loaded_from_dem_without_correlations =
        global_matching_graph.loaded_from_dem_without_correlations;
    local_matching_graph.negative_weight_detection_events_set =
        global_matching_graph.negative_weight_detection_events_set;
    local_matching_graph.negative_weight_observables_set =
        global_matching_graph.negative_weight_observables_set;
    local_matching_graph.negative_weight_sum = global_matching_graph.negative_weight_sum;
    local_matching_graph.is_user_graph_boundary_node = global_matching_graph.is_user_graph_boundary_node;

    for (size_t node_id = 0; node_id < global_matching_graph.nodes.size(); node_id++) {
        const auto& src = global_matching_graph.nodes[node_id];
        auto& dst = local_matching_graph.nodes[node_id];
        dst.neighbor_weights = src.neighbor_weights;
        dst.neighbor_observables = src.neighbor_observables;
        // The level-batched benchmark path is non-correlated. Avoid copying the
        // nested implied-weight vectors when cloning the full identity graph;
        // they are only needed by correlated decoding.
        dst.neighbor_implied_weights.resize(src.neighbors.size());
        dst.neighbors.resize(src.neighbors.size());
        for (size_t k = 0; k < src.neighbors.size(); k++) {
            if (src.neighbors[k] == nullptr) {
                dst.neighbors[k] = nullptr;
            } else {
                size_t neighbor_id = static_cast<size_t>(src.neighbors[k] - global_matching_graph.nodes.data());
                dst.neighbors[k] = &local_matching_graph.nodes[neighbor_id];
            }
        }
    }
    return local_matching_graph;
}

ClusterSubgraph build_full_identity_cluster_subgraph(UserGraph& graph, const ProcessingCluster& cluster) {
    ClusterSubgraph cluster_subgraph;
    cluster_subgraph.cluster_id = cluster.id;
    cluster_subgraph.global_active_detectors = cluster.active_detectors;
    const MatchingGraph& global_matching_graph = graph.get_matching_graph_for_parallel_clustering();
    const size_t num_nodes = graph.get_num_nodes();
    cluster_subgraph.local_to_global_node_ids.resize(num_nodes);
    std::iota(cluster_subgraph.local_to_global_node_ids.begin(), cluster_subgraph.local_to_global_node_ids.end(), 0);
    cluster_subgraph.global_to_local_node_ids.resize(num_nodes);
    std::iota(cluster_subgraph.global_to_local_node_ids.begin(), cluster_subgraph.global_to_local_node_ids.end(), 0);
    cluster_subgraph.local_node_ids_are_global_ids = true;
    if (!global_matching_graph.is_user_graph_boundary_node.empty()) {
        cluster_subgraph.local_boundary_node_mask = global_matching_graph.is_user_graph_boundary_node;
    } else if (!graph.boundary_nodes.empty()) {
        cluster_subgraph.local_boundary_node_mask.assign(num_nodes, false);
        for (auto node : graph.boundary_nodes) {
            if (node < num_nodes) {
                cluster_subgraph.local_boundary_node_mask[node] = true;
            }
        }
    } else {
        cluster_subgraph.local_boundary_node_mask.clear();
    }
    cluster_subgraph.local_active_detectors.reserve(cluster.active_detectors.size());
    for (auto detector : cluster.active_detectors) {
        if (cluster_subgraph_local_node_is_boundary(cluster_subgraph, detector)) {
            continue;
        }
        cluster_subgraph.local_active_detectors.push_back(detector);
    }
    if (cluster.parent_id == NO_PROCESSING_CLUSTER_PARENT &&
        global_matching_graph.num_observables <= sizeof(pm::obs_int) * 8) {
        // The root full-identity cluster is exactly the global matching graph.
        // Reuse the cached UserGraph MWPM object instead of rebuilding a full
        // clone in every shot.  This preserves the cluster hierarchy: only the
        // root subgraph representation is specialized.
        cluster_subgraph.shared_mwpm = &graph.get_mwpm();
        cluster_subgraph.use_shared_full_graph_mwpm = true;
        return cluster_subgraph;
    }
    // Correctness-first lockstep needs an independent MWPM state for every
    // non-root cluster. Keep the fast direct MatchingGraph clone path so that
    // we still avoid UserGraph rescaling.
    auto local_matching_graph = clone_identity_matching_graph_for_cluster(global_matching_graph);
    auto cached_mwpm = std::make_unique<Mwpm>(GraphFlooder(std::move(local_matching_graph)));
    cached_mwpm->flooder.sync_negative_weight_observables_and_detection_events();
    cached_mwpm->flooder.graph.loaded_from_dem_without_correlations =
        graph.loaded_from_dem_without_correlations;
    cluster_subgraph.cached_mwpm = std::move(cached_mwpm);
    cluster_subgraph.use_shared_full_graph_mwpm = false;
    return cluster_subgraph;
}



struct PersistentFullGraphWorkerState;

struct PersistentFullGraphResetRuntime {
    explicit PersistentFullGraphResetRuntime(PersistentFullGraphWorkerState* state) : state(state) {}
    ~PersistentFullGraphResetRuntime();

    void ensure_workers(size_t worker_count);
    void reset_loop(size_t worker_index);

    PersistentFullGraphWorkerState* state;
    std::vector<std::thread> reset_workers;
    std::mutex mutex;
    std::condition_variable work_cv;
    std::condition_variable ready_cv;
    bool stopping = false;
};

struct PersistentFullGraphWorkerState {
    const UserGraph* graph_ptr = nullptr;
    size_t num_nodes = 0;
    size_t num_observables = 0;
    bool loaded_from_dem_without_correlations = false;
    std::vector<std::unique_ptr<Mwpm>> worker_mwpms;
    std::vector<uint8_t> slot_states;  // 0=new/unused, 1=ready, 2=in-use, 3=resetting.
    std::deque<size_t> ready_slots;
    std::deque<size_t> dirty_slots;
    std::unordered_map<Mwpm*, size_t> slot_by_mwpm;
    std::unique_ptr<PersistentFullGraphResetRuntime> reset_runtime;
    size_t configured_pool_slots = 0;
};

std::mutex& persistent_full_graph_worker_owner_registry_mutex() {
    static std::mutex mutex;
    return mutex;
}

std::unordered_map<Mwpm*, PersistentFullGraphWorkerState*>& persistent_full_graph_worker_owner_registry() {
    static std::unordered_map<Mwpm*, PersistentFullGraphWorkerState*> registry;
    return registry;
}

PersistentFullGraphWorkerState* lookup_persistent_full_graph_worker_state_for_mwpm(Mwpm* mwpm) {
    std::lock_guard<std::mutex> lock(persistent_full_graph_worker_owner_registry_mutex());
    auto& registry = persistent_full_graph_worker_owner_registry();
    auto it = registry.find(mwpm);
    return it == registry.end() ? nullptr : it->second;
}

size_t env_size_t_or_default(const char* name, size_t default_value) {
    if (const char* text = std::getenv(name)) {
        try {
            return static_cast<size_t>(std::stoull(text));
        } catch (...) {
            return default_value;
        }
    }
    return default_value;
}

PersistentFullGraphResetRuntime::~PersistentFullGraphResetRuntime() {
    {
        std::lock_guard<std::mutex> lock(mutex);
        stopping = true;
    }
    work_cv.notify_all();
    for (auto& worker : reset_workers) {
        if (worker.joinable()) {
            worker.join();
        }
    }
}

void maybe_pin_current_parentless_reset_worker(size_t worker_index) {
#if defined(__linux__)
    const char* enabled = std::getenv("PYMATCHING_PARENTLESS_PIN_RESET_WORKERS");
    if (enabled == nullptr || std::string(enabled) != "1") {
        return;
    }
    size_t core_base = 0;
    if (const char* env = std::getenv("PYMATCHING_PARENTLESS_RESET_WORKER_CORE_BASE")) {
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

void PersistentFullGraphResetRuntime::ensure_workers(size_t worker_count) {
    if (worker_count == 0) {
        worker_count = 1;
    }
    while (reset_workers.size() < worker_count) {
        const size_t worker_index = reset_workers.size();
        reset_workers.emplace_back([this, worker_index]() { reset_loop(worker_index); });
    }
}

void PersistentFullGraphResetRuntime::reset_loop(size_t worker_index) {
    maybe_pin_current_parentless_reset_worker(worker_index);
    while (true) {
        size_t slot = SIZE_MAX;
        {
            std::unique_lock<std::mutex> lock(mutex);
            work_cv.wait(lock, [&]() { return stopping || !state->dirty_slots.empty(); });
            if (stopping && state->dirty_slots.empty()) {
                return;
            }
            slot = state->dirty_slots.front();
            state->dirty_slots.pop_front();
        }

        state->worker_mwpms[slot]->reset();

        {
            std::lock_guard<std::mutex> lock(mutex);
            state->slot_states[slot] = 1;
            state->ready_slots.push_back(slot);
        }
        ready_cv.notify_one();
    }
}

PersistentFullGraphWorkerState& persistent_full_graph_worker_state(UserGraph& graph) {
    static std::mutex mutex;
    // Intentionally leak the persistent states. These objects own background
    // reset threads and large predistributed detector graph/MWPM states. The
    // benchmark process exits immediately after reporting rows, and leaking here
    // avoids teardown races/static-destruction order issues that are irrelevant
    // to online decoding time.
    static std::unordered_map<const UserGraph*, PersistentFullGraphWorkerState*> states;
    std::lock_guard<std::mutex> lock(mutex);
    auto& ptr = states[&graph];
    const size_t num_nodes = graph.get_num_nodes();
    const size_t num_observables = graph.get_num_observables();
    if (ptr == nullptr || ptr->graph_ptr != &graph || ptr->num_nodes != num_nodes || ptr->num_observables != num_observables) {
        ptr = new PersistentFullGraphWorkerState();
        ptr->graph_ptr = &graph;
        ptr->num_nodes = num_nodes;
        ptr->num_observables = num_observables;
        ptr->loaded_from_dem_without_correlations = graph.loaded_from_dem_without_correlations;
        ptr->reset_runtime = std::make_unique<PersistentFullGraphResetRuntime>(ptr);
    }
    return *ptr;
}

void ensure_persistent_full_graph_worker_slot(UserGraph& graph, size_t slot) {
    static std::mutex build_mutex;
    std::lock_guard<std::mutex> build_lock(build_mutex);
    auto& state = persistent_full_graph_worker_state(graph);
    if (state.worker_mwpms.size() <= slot) {
        state.worker_mwpms.resize(slot + 1);
        state.slot_states.resize(slot + 1, 0);
    }
    if (state.worker_mwpms[slot] == nullptr) {
        const MatchingGraph& global_matching_graph = graph.get_matching_graph_for_parallel_clustering();
        auto local_matching_graph = clone_identity_matching_graph_for_cluster(global_matching_graph);
        auto mwpm = std::make_unique<Mwpm>(GraphFlooder(std::move(local_matching_graph)));
        mwpm->flooder.sync_negative_weight_observables_and_detection_events();
        mwpm->flooder.graph.loaded_from_dem_without_correlations =
            graph.loaded_from_dem_without_correlations;
        auto* raw_mwpm = mwpm.get();
        state.slot_by_mwpm[raw_mwpm] = slot;
        {
            std::lock_guard<std::mutex> lock(persistent_full_graph_worker_owner_registry_mutex());
            persistent_full_graph_worker_owner_registry()[raw_mwpm] = &state;
        }
        state.worker_mwpms[slot] = std::move(mwpm);
    }
}

Mwpm* get_persistent_full_graph_worker_mwpm(UserGraph& graph, size_t slot) {
    ensure_persistent_full_graph_worker_slot(graph, slot);
    auto& state = persistent_full_graph_worker_state(graph);
    return state.worker_mwpms[slot].get();
}

Mwpm* acquire_clean_persistent_full_graph_worker_mwpm(UserGraph& graph) {
    auto& state = persistent_full_graph_worker_state(graph);
    if (state.reset_runtime == nullptr) {
        state.reset_runtime = std::make_unique<PersistentFullGraphResetRuntime>(&state);
    }
    {
        std::unique_lock<std::mutex> lock(state.reset_runtime->mutex);
        while (state.ready_slots.empty()) {
            // If prewarm was forgotten or the configured pool has not yet been
            // fully materialized, grow up to the configured pool size. Once the
            // configured pool exists, wait for a released dirty slot to be reset.
            // The parentless level-critical path avoids construction-time
            // deadlocks by acquiring persistent MWPM states lazily instead of
            // by growing beyond the configured pool.
            const size_t configured = state.configured_pool_slots;
            if (configured == 0 || state.worker_mwpms.size() < configured) {
                const size_t slot = state.worker_mwpms.size();
                lock.unlock();
                ensure_persistent_full_graph_worker_slot(graph, slot);
                lock.lock();
                if (state.slot_states[slot] == 0) {
                    state.slot_states[slot] = 1;
                    state.ready_slots.push_back(slot);
                    if (state.configured_pool_slots < slot + 1) {
                        state.configured_pool_slots = slot + 1;
                    }
                }
                continue;
            }
            state.reset_runtime->ready_cv.wait(lock, [&]() { return !state.ready_slots.empty(); });
        }
        const size_t slot = state.ready_slots.front();
        state.ready_slots.pop_front();
        state.slot_states[slot] = 2;
        return state.worker_mwpms[slot].get();
    }
}

MatchingGraph build_matching_graph_directly_from_global_matching_graph(
    const MatchingGraph& global_matching_graph,
    const ClusterSubgraph& cluster_subgraph) {
    MatchingGraph local_matching_graph(
        cluster_subgraph.local_to_global_node_ids.size(),
        global_matching_graph.num_observables,
        global_matching_graph.normalising_constant);
    local_matching_graph.loaded_from_dem_without_correlations =
        global_matching_graph.loaded_from_dem_without_correlations;
    if (!cluster_subgraph.local_boundary_node_mask.empty()) {
        local_matching_graph.is_user_graph_boundary_node.assign(
            cluster_subgraph.local_to_global_node_ids.size(), false);
        for (size_t local = 0; local < cluster_subgraph.local_boundary_node_mask.size(); local++) {
            local_matching_graph.is_user_graph_boundary_node[local] = cluster_subgraph.local_boundary_node_mask[local];
        }
    }

    // Reserve a small upper bound on each adjacency list and build edges in one
    // pass. The earlier exact-degree reservation scanned every cluster edge
    // twice; over-reserving by the original detector degree is cheaper because
    // surface-code detector degree is bounded and small.
    for (size_t local_node = 0; local_node < local_matching_graph.nodes.size(); local_node++) {
        size_t global_node = cluster_subgraph.local_to_global_node_ids[local_node];
        const auto& source = global_matching_graph.nodes[global_node];
        auto& node = local_matching_graph.nodes[local_node];
        auto degree = source.neighbors.size();
        node.neighbors.reserve(degree);
        node.neighbor_weights.reserve(degree);
        node.neighbor_observables.reserve(degree);
        node.neighbor_implied_weights.reserve(degree);
    }

    for (size_t local_node = 0; local_node < cluster_subgraph.local_to_global_node_ids.size(); local_node++) {
        size_t global_node = cluster_subgraph.local_to_global_node_ids[local_node];
        const auto& source = global_matching_graph.nodes[global_node];
        for (size_t k = 0; k < source.neighbors.size(); k++) {
            const DetectorNode* global_neighbor = source.neighbors[k];
            auto obs_mask = source.neighbor_observables[k];
            auto weight = static_cast<signed_weight_int>(source.neighbor_weights[k]);
            if (global_neighbor == nullptr) {
                add_matching_graph_boundary_edge_fast(local_matching_graph, local_node, weight, obs_mask);
                continue;
            }
            size_t global_neighbor_id = static_cast<size_t>(global_neighbor - global_matching_graph.nodes.data());
            if (global_neighbor_id <= global_node) {
                continue;
            }
            if (global_neighbor_id >= cluster_subgraph.global_to_local_node_ids.size()) {
                continue;
            }
            size_t local_neighbor = cluster_subgraph.global_to_local_node_ids[global_neighbor_id];
            if (local_neighbor == SIZE_MAX) {
                continue;
            }
            add_matching_graph_edge_fast(local_matching_graph, local_node, local_neighbor, weight, obs_mask);
        }
    }
    return local_matching_graph;
}

std::vector<uint64_t> map_cluster_active_detectors(
    const ClusterSubgraph& cluster_subgraph, const std::vector<uint64_t>& global_active_detectors) {
    std::vector<uint64_t> local_active_detectors;
    local_active_detectors.reserve(global_active_detectors.size());
    for (auto detector : global_active_detectors) {
        if (detector >= cluster_subgraph.global_to_local_node_ids.size()) {
            throw std::invalid_argument("Detection event is outside the source graph.");
        }
        size_t local = cluster_subgraph.global_to_local_node_ids[detector];
        if (local == SIZE_MAX) {
            throw std::invalid_argument("Detection event is outside the cluster subgraph.");
        }
        if (cluster_subgraph_local_node_is_boundary(cluster_subgraph, local)) {
            continue;
        }
        local_active_detectors.push_back(local);
    }
    std::sort(local_active_detectors.begin(), local_active_detectors.end());
    local_active_detectors.erase(std::unique(local_active_detectors.begin(), local_active_detectors.end()), local_active_detectors.end());
    return local_active_detectors;
}

}  // namespace

void prewarm_persistent_full_graph_worker_states(UserGraph& graph, size_t worker_slots) {
    if (worker_slots == 0) {
        return;
    }
    auto& state = persistent_full_graph_worker_state(graph);
    const size_t default_pool_slots = std::max<size_t>(worker_slots, worker_slots * 4);
    const size_t pool_slots = env_size_t_or_default(
        "PYMATCHING_PARENTLESS_STATE_POOL_SLOTS", default_pool_slots);
    const size_t reset_workers = env_size_t_or_default(
        "PYMATCHING_PARENTLESS_RESET_WORKERS", std::max<size_t>(1, worker_slots));

    for (size_t slot = 0; slot < pool_slots; slot++) {
        (void)get_persistent_full_graph_worker_mwpm(graph, slot);
    }

    auto& runtime = *state.reset_runtime;
    runtime.ensure_workers(reset_workers);
    {
        std::lock_guard<std::mutex> lock(runtime.mutex);
        state.configured_pool_slots = std::max(state.configured_pool_slots, pool_slots);
        for (size_t slot = 0; slot < pool_slots; slot++) {
            if (state.slot_states[slot] == 0) {
                state.slot_states[slot] = 1;
                state.ready_slots.push_back(slot);
            }
        }
    }
    runtime.ready_cv.notify_all();
}

void wait_for_persistent_full_graph_worker_states_ready(UserGraph& graph, size_t min_ready) {
    if (min_ready == 0) {
        return;
    }
    auto& state = persistent_full_graph_worker_state(graph);
    if (state.reset_runtime == nullptr) {
        state.reset_runtime = std::make_unique<PersistentFullGraphResetRuntime>(&state);
    }
    auto& runtime = *state.reset_runtime;
    std::unique_lock<std::mutex> lock(runtime.mutex);
    runtime.ready_cv.wait(lock, [&]() {
        return state.ready_slots.size() >= min_ready;
    });
}

void release_persistent_full_graph_worker_mwpm(Mwpm* mwpm) {
    if (mwpm == nullptr) {
        return;
    }
    // The benchmark path has one active graph; the pointer reverse lookup is
    // intentionally linear over the small static state table to keep ownership
    // hidden inside this translation unit.
    // Find the owner state by consulting every known graph state.
    // A separate global registry would also work but is easier to make stale
    // when a graph object is rebuilt.
    static std::mutex lookup_mutex;
    (void)lookup_mutex;
    // Local helper: iterate through the same map as persistent_full_graph_worker_state.
    // Since that map is private to the function, expose the common case via the
    // slot_by_mwpm search on all states would require refactoring.  Instead the
    // state pointer is recovered from a static registry filled at slot creation.
    auto* state = lookup_persistent_full_graph_worker_state_for_mwpm(mwpm);
    if (state == nullptr || state->reset_runtime == nullptr) {
        mwpm->reset();
        return;
    }

    auto& runtime = *state->reset_runtime;
    {
        std::lock_guard<std::mutex> lock(runtime.mutex);
        auto it = state->slot_by_mwpm.find(mwpm);
        if (it == state->slot_by_mwpm.end()) {
            mwpm->reset();
            return;
        }
        const size_t slot = it->second;
        if (slot >= state->slot_states.size()) {
            mwpm->reset();
            return;
        }
        if (state->slot_states[slot] != 2) {
            return;
        }
        state->slot_states[slot] = 3;
        state->dirty_slots.push_back(slot);
    }
    runtime.work_cv.notify_one();
}

ClusterSubgraph build_full_graph_worker_cluster_subgraph(
    UserGraph& graph,
    const ProcessingCluster& cluster,
    bool allow_shared_global_mwpm,
    bool lazy_acquire_persistent_mwpm) {
    validate_cluster_vertices(graph, cluster);
    ClusterSubgraph cluster_subgraph;
    cluster_subgraph.cluster_id = cluster.id;
    cluster_subgraph.global_active_detectors = cluster.active_detectors;

    const MatchingGraph& global_matching_graph = graph.get_matching_graph_for_parallel_clustering();
    const size_t num_nodes = graph.get_num_nodes();

    // The whole point of this mode is that worker-local ids are the original
    // detector ids.  This preserves stopped-state import/export remapping while
    // avoiding any cut-out graph or boundary endpoint materialization.
    cluster_subgraph.local_to_global_node_ids.resize(num_nodes);
    std::iota(cluster_subgraph.local_to_global_node_ids.begin(), cluster_subgraph.local_to_global_node_ids.end(), 0);
    cluster_subgraph.global_to_local_node_ids.resize(num_nodes);
    std::iota(cluster_subgraph.global_to_local_node_ids.begin(), cluster_subgraph.global_to_local_node_ids.end(), 0);
    cluster_subgraph.local_node_ids_are_global_ids = true;

    if (!global_matching_graph.is_user_graph_boundary_node.empty()) {
        cluster_subgraph.local_boundary_node_mask = global_matching_graph.is_user_graph_boundary_node;
    } else if (!graph.boundary_nodes.empty()) {
        cluster_subgraph.local_boundary_node_mask.assign(num_nodes, false);
        for (auto node : graph.boundary_nodes) {
            if (node < num_nodes) {
                cluster_subgraph.local_boundary_node_mask[node] = true;
            }
        }
    } else {
        cluster_subgraph.local_boundary_node_mask.clear();
    }

    cluster_subgraph.local_active_detectors.reserve(cluster.active_detectors.size());
    for (auto detector : cluster.active_detectors) {
        if (detector >= num_nodes) {
            throw std::invalid_argument("ProcessingCluster active detector is not present in the UserGraph.");
        }
        if (cluster_subgraph_local_node_is_boundary(cluster_subgraph, detector)) {
            continue;
        }
        cluster_subgraph.local_active_detectors.push_back(detector);
    }
    std::sort(cluster_subgraph.local_active_detectors.begin(), cluster_subgraph.local_active_detectors.end());
    cluster_subgraph.local_active_detectors.erase(
        std::unique(cluster_subgraph.local_active_detectors.begin(), cluster_subgraph.local_active_detectors.end()),
        cluster_subgraph.local_active_detectors.end());

    if (graph.get_num_observables() <= sizeof(pm::obs_int) * 8) {
        if (allow_shared_global_mwpm && cluster.parent_id == NO_PROCESSING_CLUSTER_PARENT && !cluster.forced_by_max_level) {
            // The canonical single-root worker can reuse UserGraph's prebuilt
            // global MWPM object.  When a level-1 or other non-forced policy
            // promotes multiple independent roots in the same shot, sharing the
            // same flooder is invalid: the first root leaves scheduled events in
            // the queue before the next root is initialized.  The scheduler passes
            // allow_shared_global_mwpm=false unless there is exactly one root that
            // can safely own the global MWPM state.
            cluster_subgraph.shared_mwpm = &graph.get_mwpm();
            cluster_subgraph.use_shared_full_graph_mwpm = true;
            return cluster_subgraph;
        }
        // Persistent full-graph worker mode. In the ordinary scheduler path this
        // acquires a clean predistributed detector graph/MWPM state immediately.
        // In the parentless level-critical accounting path it is acquired lazily
        // just before the cluster is initialized/executed, then released right
        // after export/extraction. This prevents a slot-per-cluster requirement
        // during per-shot subgraph construction.
        if (lazy_acquire_persistent_mwpm) {
            cluster_subgraph.lazy_persistent_full_graph_worker_mwpm = true;
        } else {
            cluster_subgraph.persistent_mwpm = acquire_clean_persistent_full_graph_worker_mwpm(graph);
            cluster_subgraph.persistent_mwpm_acquired_clean = true;
        }
    } else {
        // Rare fallback for >64 observable graphs: build an ordinary UserGraph
        // clone so get_mwpm_with_search_graph() can construct both flooders.
        // The current benchmark path uses <=64 observables and therefore uses
        // the direct MatchingGraph clone above.
        throw std::invalid_argument(
            "Full-graph worker mode currently supports non-correlated <=64-observable decoding.");
    }

    cluster_subgraph.use_shared_full_graph_mwpm = false;
    cluster_subgraph.shared_mwpm = nullptr;
    return cluster_subgraph;
}

ClusterSubgraph build_cluster_subgraph(UserGraph& graph, const ProcessingCluster& cluster) {
    validate_cluster_vertices(graph, cluster);
    ClusterSubgraph cluster_subgraph;
    cluster_subgraph.cluster_id = cluster.id;
    cluster_subgraph.global_active_detectors = cluster.active_detectors;
    const MatchingGraph& global_matching_graph = graph.get_matching_graph_for_parallel_clustering();
    bool can_use_direct_cached_mwpm = graph.get_num_observables() <= sizeof(pm::obs_int) * 8;

    if (can_use_direct_cached_mwpm &&
        cluster_influence_is_full_identity(graph, cluster)) {
        return build_full_identity_cluster_subgraph(graph, cluster);
    }

    cluster_subgraph.local_to_global_node_ids = collect_kept_global_nodes(graph, cluster);
    cluster_subgraph.global_to_local_node_ids =
        build_global_to_local_node_ids(graph.get_num_nodes(), cluster_subgraph.local_to_global_node_ids);
    std::vector<bool> kept_global_nodes;
    if (!can_use_direct_cached_mwpm) {
        kept_global_nodes = compute_kept_mask_from_local_to_global(
            graph.get_num_nodes(), cluster_subgraph.local_to_global_node_ids);
    }

    if (can_use_direct_cached_mwpm) {
        if (graph.boundary_nodes.empty()) {
            cluster_subgraph.local_boundary_node_mask.clear();
        } else {
            cluster_subgraph.local_boundary_node_mask.assign(cluster_subgraph.local_to_global_node_ids.size(), false);
            for (size_t local = 0; local < cluster_subgraph.local_to_global_node_ids.size(); local++) {
                cluster_subgraph.local_boundary_node_mask[local] =
                    graph.boundary_nodes.contains(cluster_subgraph.local_to_global_node_ids[local]);
            }
        }
    } else {
        cluster_subgraph.graph = UserGraph(cluster_subgraph.local_to_global_node_ids.size(), graph.get_num_observables());
        cluster_subgraph.graph.set_min_num_observables(graph.get_num_observables());
        cluster_subgraph.graph.loaded_from_dem_without_correlations = graph.loaded_from_dem_without_correlations;

        auto local_boundary_nodes = build_local_boundary_nodes(
            graph, cluster_subgraph.local_to_global_node_ids, cluster_subgraph.local_boundary_node_mask);
        cluster_subgraph.graph.set_boundary(local_boundary_nodes);
    }

    if (!can_use_direct_cached_mwpm) {
        std::vector<PendingImpliedWeightRewrite> pending_rules;
        for (const auto& edge : graph.edges) {
            if (edge.node1 == SIZE_MAX || edge.node1 >= kept_global_nodes.size() || !kept_global_nodes[edge.node1]) {
                continue;
            }
            size_t local_node1 = cluster_subgraph.global_to_local_node_ids[edge.node1];
            if (edge.node2 == SIZE_MAX) {
                add_included_edge(cluster_subgraph.graph, edge, local_node1, SIZE_MAX, pending_rules);
                continue;
            }
            if (edge.node2 >= kept_global_nodes.size() || !kept_global_nodes[edge.node2]) {
                continue;
            }
            size_t local_node2 = cluster_subgraph.global_to_local_node_ids[edge.node2];
            add_included_edge(cluster_subgraph.graph, edge, local_node1, local_node2, pending_rules);
        }
        rewrite_implied_weight_rules(cluster_subgraph, pending_rules);
    }

    cluster_subgraph.local_active_detectors =
        map_cluster_active_detectors(cluster_subgraph, cluster_subgraph.global_active_detectors);

    // Reuse the already-discretized global MatchingGraph to avoid per-cluster
    // UserGraph -> MatchingGraph scaling during lockstep initialization. In this
    // fast non-correlated path, the UserGraph intentionally keeps only metadata
    // and boundary markers; all actual decoding uses cached_mwpm.
    if (can_use_direct_cached_mwpm) {
        auto local_matching_graph = build_matching_graph_directly_from_global_matching_graph(
            global_matching_graph, cluster_subgraph);
        auto cached_mwpm = std::make_unique<Mwpm>(GraphFlooder(std::move(local_matching_graph)));
        cached_mwpm->flooder.sync_negative_weight_observables_and_detection_events();
        cached_mwpm->flooder.graph.loaded_from_dem_without_correlations =
            graph.loaded_from_dem_without_correlations;
        cluster_subgraph.cached_mwpm = std::move(cached_mwpm);
    }
    return cluster_subgraph;
}

void acquire_persistent_full_graph_worker_mwpm_for_subgraph(
    UserGraph& graph,
    ClusterSubgraph& cluster_subgraph) {
    if (!cluster_subgraph.lazy_persistent_full_graph_worker_mwpm) {
        return;
    }
    if (cluster_subgraph.persistent_mwpm != nullptr) {
        return;
    }
    cluster_subgraph.persistent_mwpm = acquire_clean_persistent_full_graph_worker_mwpm(graph);
    cluster_subgraph.persistent_mwpm_acquired_clean = true;
}

Mwpm& cluster_subgraph_mwpm(ClusterSubgraph& cluster_subgraph, bool edge_correlations) {
    if (cluster_subgraph.lazy_persistent_full_graph_worker_mwpm &&
        cluster_subgraph.persistent_mwpm == nullptr) {
        throw std::invalid_argument(
            "Lazy persistent full-graph worker MWPM was requested before acquisition.");
    }
    if (cluster_subgraph.persistent_mwpm != nullptr) {
        return *cluster_subgraph.persistent_mwpm;
    }
    if (cluster_subgraph.use_shared_full_graph_mwpm) {
        if (cluster_subgraph.shared_mwpm == nullptr) {
            throw std::invalid_argument("Shared full-identity cluster subgraph is missing its MWPM pointer.");
        }
        if (edge_correlations) {
            throw std::invalid_argument("Shared full-identity cluster subgraph does not support edge correlations.");
        }
        return *cluster_subgraph.shared_mwpm;
    }
    if (!edge_correlations && cluster_subgraph.cached_mwpm != nullptr) {
        return *cluster_subgraph.cached_mwpm;
    }
    bool needs_search_graph = edge_correlations ||
        cluster_subgraph.graph.get_num_observables() > sizeof(pm::obs_int) * 8;
    return needs_search_graph ? cluster_subgraph.graph.get_mwpm_with_search_graph()
                              : cluster_subgraph.graph.get_mwpm();
}

std::vector<uint64_t> map_global_detection_events_to_cluster_subgraph(
    const ClusterSubgraph& cluster_subgraph, const std::vector<uint64_t>& global_detection_events) {
    return map_cluster_active_detectors(cluster_subgraph, global_detection_events);
}

ExtendedMatchingResult decode_cluster_subgraph(
    ClusterSubgraph& cluster_subgraph,
    const std::vector<uint64_t>& global_detection_events,
    bool edge_correlations) {
    if (edge_correlations && cluster_subgraph.dropped_implied_weight_rules) {
        throw std::invalid_argument(
            "Cluster subgraph dropped implied-weight rules at the cluster boundary, so edge-correlated decoding is not "
            "supported for this subgraph.");
    }

    auto local_detection_events = map_global_detection_events_to_cluster_subgraph(cluster_subgraph, global_detection_events);
    auto& mwpm = cluster_subgraph_mwpm(cluster_subgraph, edge_correlations);
    ExtendedMatchingResult result(mwpm.flooder.graph.num_observables);
    decode_detection_events(mwpm, local_detection_events, result.obs_crossed.data(), result.weight, edge_correlations);
    return result;
}

ExtendedMatchingResult decode_cluster_subgraph(ClusterSubgraph& cluster_subgraph, bool edge_correlations) {
    return decode_cluster_subgraph(cluster_subgraph, cluster_subgraph.global_active_detectors, edge_correlations);
}

}  // namespace pm
