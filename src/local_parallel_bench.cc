#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>

#include "pymatching/sparse_blossom/driver/mwpm_decoding.h"
#include "pymatching/sparse_blossom/driver/user_graph.h"
#include "pymatching/sparse_blossom/parallel/lockstep_scheduler.h"
#include "pymatching/sparse_blossom/parallel/growing_only_clusterer.h"

using Clock = std::chrono::steady_clock;

struct Shot { uint64_t obs_mask; std::vector<uint64_t> hits; };
struct CaseData { pm::UserGraph graph; std::vector<Shot> shots; size_t rounds; size_t distance; double noise; size_t num_edges; };

bool file_has_content(const std::string& path) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    return in.good() && in.tellg() > 0;
}


bool env_flag_enabled(const char* name);

bool eventcount_wall_timing_disabled() {
    return env_flag_enabled("PYMATCHING_DISABLE_WALL_TIMING") ||
           env_flag_enabled("PYMATCHING_EVENTCOUNT_ONLY_NO_WALL_TIMING");
}

Clock::time_point eventcount_optional_now() {
    return eventcount_wall_timing_disabled() ? Clock::time_point{} : Clock::now();
}


template <typename T>
uint64_t value_or_zero(const std::vector<T>& values, size_t index) {
    return index < values.size() ? static_cast<uint64_t>(values[index]) : 0;
}

void append_cluster_event_csv_rows(
    const std::string& output_path,
    const std::string& case_name,
    const CaseData& data,
    int rep_index,
    size_t shot_index,
    const pm::LockstepDecodeResult& result) {
    const bool write_header = !file_has_content(output_path);
    std::ofstream out(output_path, std::ios::app);
    if (!out) {
        throw std::runtime_error("failed to open cluster event csv output: " + output_path);
    }
    if (write_header) {
        out << "case,distance,rounds,noise,rep,shot,cluster_id,level,active_detectors,"
               "influence_vertices,boundary_vertices,diameter_bound,buffer_bound,diameter,"
               "nearest_boundary_match_distance,nearest_external_detector_distance,forced_by_max_level,"
               "event_count,critical_path_event_count,has_cutoff_algorithmic_time,cutoff_algorithmic_time,"
               "prior_level_critical_event_count,overlap_event_base,"
               "events_at_or_before_cutoff,events_after_cutoff,stop_algorithmic_time,critical_path_algorithmic_time\n";
    }
    for (size_t cluster_id = 0; cluster_id < result.hierarchy.clusters.size(); cluster_id++) {
        const auto& cluster = result.hierarchy.clusters[cluster_id];
        out << case_name
            << "," << data.distance
            << "," << data.rounds
            << "," << data.noise
            << "," << rep_index
            << "," << shot_index
            << "," << cluster.id
            << "," << cluster.level
            << "," << cluster.active_detectors.size()
            << "," << cluster.influence_vertices.size()
            << "," << cluster.boundary_endpoint_vertices.size()
            << "," << cluster.diameter_bound
            << "," << cluster.buffer_bound
            << "," << cluster.diameter
            << "," << cluster.nearest_boundary_match_distance
            << "," << cluster.nearest_external_detector_distance
            << "," << (cluster.forced_by_max_level ? 1 : 0)
            << "," << value_or_zero(result.cluster_step_event_counts, cluster_id)
            << "," << value_or_zero(result.cluster_critical_path_event_counts, cluster_id)
            << "," << value_or_zero(result.cluster_has_cutoff_algorithmic_time, cluster_id)
            << "," << value_or_zero(result.cluster_cutoff_algorithmic_times, cluster_id)
            << "," << value_or_zero(result.cluster_prior_level_critical_event_counts, cluster_id)
            << "," << value_or_zero(result.cluster_overlap_event_bases, cluster_id)
            << "," << value_or_zero(result.cluster_events_at_or_before_cutoff, cluster_id)
            << "," << value_or_zero(result.cluster_events_after_cutoff, cluster_id)
            << "," << value_or_zero(result.cluster_stop_algorithmic_times, cluster_id)
            << "," << value_or_zero(result.cluster_critical_path_algorithmic_times, cluster_id)
            << "\n";
    }

}

std::string method_label_from_env(const char* fallback) {
    if (const char* e = std::getenv("PYMATCHING_EVENTCOUNT_METHOD_LABEL")) {
        if (*e != '\0') {
            return std::string(e);
        }
    }
    return std::string(fallback);
}

std::string join_u64_semicolon(std::vector<uint64_t> values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    std::ostringstream out;
    for (size_t k = 0; k < values.size(); k++) {
        if (k) out << ";";
        out << values[k];
    }
    return out.str();
}

struct ClusterDetailCsvRow {
    size_t cluster_id = 0;
    size_t level = 0;
    std::vector<uint64_t> active_detectors;
    pm::cumulative_time_int diameter_bound = 0;
    pm::cumulative_time_int buffer_bound = 0;
    pm::cumulative_time_int diameter = 0;
    pm::cumulative_time_int nearest_boundary_match_distance = 0;
    pm::cumulative_time_int nearest_external_detector_distance = 0;
    bool forced_by_max_level = false;
    uint64_t event_count = 0;
    uint64_t critical_path_event_count = 0;
    uint64_t has_cutoff_algorithmic_time = 0;
    uint64_t cutoff_algorithmic_time = 0;
    uint64_t prior_level_critical_event_count = 0;
    uint64_t overlap_event_base = 0;
    uint64_t events_at_or_before_cutoff = 0;
    uint64_t events_after_cutoff = 0;
    uint64_t stop_algorithmic_time = 0;
    uint64_t critical_path_algorithmic_time = 0;
};

struct LevelDetailSummary {
    pm::cumulative_time_int diameter_bound = 0;
    pm::cumulative_time_int buffer_bound = 0;
    bool have_bounds = false;
    size_t cluster_count = 0;
    size_t assigned_active_detector_count = 0;
    size_t max_cluster_active_detectors = 0;
    uint64_t max_cluster_event_count = 0;
    uint64_t max_cluster_critical_path_event_count = 0;
};

void append_cluster_detail_csv_rows(
    const std::string& members_path,
    const std::string& assignments_path,
    const std::string& level_params_path,
    const std::string& method,
    const std::string& case_name,
    const CaseData& data,
    int rep_index,
    size_t shot_index,
    const std::vector<ClusterDetailCsvRow>& rows) {
    if (!members_path.empty()) {
        const bool write_header = !file_has_content(members_path);
        std::ofstream out(members_path, std::ios::app);
        if (!out) {
            throw std::runtime_error("failed to open cluster members csv output: " + members_path);
        }
        if (write_header) {
            out << "case,distance,rounds,noise,rep,shot,method,cluster_id,level,"
                   "active_detector_count,active_detectors,diameter_bound,buffer_bound,diameter,"
                   "nearest_boundary_match_distance,nearest_external_detector_distance,forced_by_max_level,"
                   "event_count,critical_path_event_count,has_cutoff_algorithmic_time,cutoff_algorithmic_time,"
                   "prior_level_critical_event_count,overlap_event_base,"
                   "events_at_or_before_cutoff,events_after_cutoff,stop_algorithmic_time,critical_path_algorithmic_time\n";
        }
        for (const auto& row : rows) {
            out << case_name
                << "," << data.distance
                << "," << data.rounds
                << "," << data.noise
                << "," << rep_index
                << "," << shot_index
                << "," << method
                << "," << row.cluster_id
                << "," << row.level
                << "," << row.active_detectors.size()
                << "," << join_u64_semicolon(row.active_detectors)
                << "," << row.diameter_bound
                << "," << row.buffer_bound
                << "," << row.diameter
                << "," << row.nearest_boundary_match_distance
                << "," << row.nearest_external_detector_distance
                << "," << (row.forced_by_max_level ? 1 : 0)
                << "," << row.event_count
                << "," << row.critical_path_event_count
                << "," << row.has_cutoff_algorithmic_time
                << "," << row.cutoff_algorithmic_time
                << "," << row.prior_level_critical_event_count
                << "," << row.overlap_event_base
                << "," << row.events_at_or_before_cutoff
                << "," << row.events_after_cutoff
                << "," << row.stop_algorithmic_time
                << "," << row.critical_path_algorithmic_time
                << "\n";
        }
    }

    if (!assignments_path.empty()) {
        const bool write_header = !file_has_content(assignments_path);
        std::ofstream out(assignments_path, std::ios::app);
        if (!out) {
            throw std::runtime_error("failed to open detector assignments csv output: " + assignments_path);
        }
        if (write_header) {
            out << "case,distance,rounds,noise,rep,shot,method,active_detector,cluster_id,level,"
                   "cluster_active_detector_count,cluster_active_detectors,diameter_bound,buffer_bound,diameter,"
                   "event_count,critical_path_event_count,has_cutoff_algorithmic_time,cutoff_algorithmic_time,"
                   "prior_level_critical_event_count,overlap_event_base,"
                   "events_at_or_before_cutoff,events_after_cutoff,stop_algorithmic_time,critical_path_algorithmic_time\n";
        }
        for (const auto& row : rows) {
            const auto cluster_active_detectors = join_u64_semicolon(row.active_detectors);
            for (auto detector : row.active_detectors) {
                out << case_name
                    << "," << data.distance
                    << "," << data.rounds
                    << "," << data.noise
                    << "," << rep_index
                    << "," << shot_index
                    << "," << method
                    << "," << detector
                    << "," << row.cluster_id
                    << "," << row.level
                    << "," << row.active_detectors.size()
                    << "," << cluster_active_detectors
                    << "," << row.diameter_bound
                    << "," << row.buffer_bound
                    << "," << row.diameter
                    << "," << row.event_count
                    << "," << row.critical_path_event_count
                    << "," << row.has_cutoff_algorithmic_time
                    << "," << row.cutoff_algorithmic_time
                    << "," << row.prior_level_critical_event_count
                    << "," << row.overlap_event_base
                    << "," << row.events_at_or_before_cutoff
                    << "," << row.events_after_cutoff
                    << "," << row.stop_algorithmic_time
                    << "," << row.critical_path_algorithmic_time
                    << "\n";
            }
        }
    }

    if (!level_params_path.empty()) {
        std::map<size_t, LevelDetailSummary> by_level;
        for (const auto& row : rows) {
            auto& summary = by_level[row.level];
            if (!summary.have_bounds) {
                summary.diameter_bound = row.diameter_bound;
                summary.buffer_bound = row.buffer_bound;
                summary.have_bounds = true;
            } else {
                summary.diameter_bound = std::max(summary.diameter_bound, row.diameter_bound);
                summary.buffer_bound = std::max(summary.buffer_bound, row.buffer_bound);
            }
            summary.cluster_count++;
            summary.assigned_active_detector_count += row.active_detectors.size();
            summary.max_cluster_active_detectors = std::max(summary.max_cluster_active_detectors, row.active_detectors.size());
            summary.max_cluster_event_count = std::max(summary.max_cluster_event_count, row.event_count);
            summary.max_cluster_critical_path_event_count = std::max(summary.max_cluster_critical_path_event_count, row.critical_path_event_count);
        }
        const bool write_header = !file_has_content(level_params_path);
        std::ofstream out(level_params_path, std::ios::app);
        if (!out) {
            throw std::runtime_error("failed to open level params csv output: " + level_params_path);
        }
        if (write_header) {
            out << "case,distance,rounds,noise,rep,shot,method,level,diameter_bound,buffer_bound,"
                   "cluster_count,assigned_active_detector_count,max_cluster_active_detectors,"
                   "max_cluster_event_count,max_cluster_critical_path_event_count\n";
        }
        for (const auto& kv : by_level) {
            const auto& summary = kv.second;
            out << case_name
                << "," << data.distance
                << "," << data.rounds
                << "," << data.noise
                << "," << rep_index
                << "," << shot_index
                << "," << method
                << "," << kv.first
                << "," << summary.diameter_bound
                << "," << summary.buffer_bound
                << "," << summary.cluster_count
                << "," << summary.assigned_active_detector_count
                << "," << summary.max_cluster_active_detectors
                << "," << summary.max_cluster_event_count
                << "," << summary.max_cluster_critical_path_event_count
                << "\n";
        }
    }
}



std::string lockstep_state_name(pm::LockstepClusterExecutionState state) {
    switch (state) {
        case pm::LockstepClusterExecutionState::RUNNING: return "RUNNING";
        case pm::LockstepClusterExecutionState::QUIESCENT_STOPPED: return "QUIESCENT_STOPPED";
        case pm::LockstepClusterExecutionState::DRAINED_WITH_LIVE_ALT_TREE: return "DRAINED_WITH_LIVE_ALT_TREE";
        case pm::LockstepClusterExecutionState::IMPORTED_INTO_PARENT: return "IMPORTED_INTO_PARENT";
    }
    return "UNKNOWN";
}

std::string join_obs_bits_from_mask(pm::obs_int mask) {
    std::vector<uint64_t> bits;
    for (size_t k = 0; k < sizeof(pm::obs_int) * 8; k++) {
        if ((mask >> k) & 1) bits.push_back(k);
    }
    return join_u64_semicolon(bits);
}

template <typename BoolLikeVector>
std::string join_bool_obs_bits(const BoolLikeVector& obs_crossed) {
    std::vector<uint64_t> bits;
    for (size_t k = 0; k < obs_crossed.size(); k++) {
        if (obs_crossed[k]) bits.push_back(k);
    }
    return join_u64_semicolon(bits);
}

std::string stopped_edge_summary(const pm::StoppedCompressedEdge& edge) {
    if (edge.loc_from == pm::NO_STOPPED_NODE && edge.loc_to == pm::NO_STOPPED_NODE) {
        return "none";
    }
    return std::to_string(edge.loc_from) + "-" + std::to_string(edge.loc_to) + "@obs=" +
           std::to_string((uint64_t)edge.obs_mask);
}

std::string stopped_region_match_summary(const pm::MwpmStoppedState& state) {
    std::vector<std::string> parts;
    for (const auto& region : state.region_states) {
        const bool has_match =
            region.match.region_id != pm::NO_STOPPED_REGION ||
            region.match.edge.loc_from != pm::NO_STOPPED_NODE ||
            region.match.edge.loc_to != pm::NO_STOPPED_NODE;
        if (!has_match) continue;
        parts.push_back(
            std::to_string(region.region_id) + "->" +
            (region.match.region_id == pm::NO_STOPPED_REGION ? std::string("none") : std::to_string(region.match.region_id)) +
            ":" + stopped_edge_summary(region.match.edge));
    }
    std::ostringstream out;
    for (size_t k = 0; k < parts.size(); k++) {
        if (k) out << ";";
        out << parts[k];
    }
    return out.str();
}

std::string stopped_detector_ownership_summary(const pm::MwpmStoppedState& state) {
    std::vector<std::string> parts;
    for (const auto& det : state.detector_states) {
        if (det.node_id == pm::NO_STOPPED_NODE) continue;
        parts.push_back(
            std::to_string(det.node_id) + "->r" +
            (det.region_id == pm::NO_STOPPED_REGION ? std::string("none") : std::to_string(det.region_id)) +
            ":src" +
            (det.reached_from_source_id == pm::NO_STOPPED_NODE ? std::string("none") : std::to_string(det.reached_from_source_id)) +
            ":obs" + std::to_string((uint64_t)det.observables_crossed_from_source));
    }
    std::ostringstream out;
    for (size_t k = 0; k < parts.size(); k++) {
        if (k) out << ";";
        out << parts[k];
    }
    return out.str();
}

size_t stopped_region_match_count(const pm::MwpmStoppedState& state) {
    size_t count = 0;
    for (const auto& region : state.region_states) {
        if (region.match.region_id != pm::NO_STOPPED_REGION ||
            region.match.edge.loc_from != pm::NO_STOPPED_NODE ||
            region.match.edge.loc_to != pm::NO_STOPPED_NODE) {
            count++;
        }
    }
    return count;
}

size_t stopped_top_region_count(const pm::MwpmStoppedState& state) {
    size_t count = 0;
    for (const auto& region : state.region_states) {
        if (region.blossom_parent_id == pm::NO_STOPPED_REGION) count++;
    }
    return count;
}

size_t stopped_blossom_child_edge_count(const pm::MwpmStoppedState& state) {
    size_t count = 0;
    for (const auto& region : state.region_states) count += region.blossom_children.size();
    return count;
}

std::vector<ClusterDetailCsvRow> make_cluster_detail_rows_from_result(const pm::LockstepDecodeResult& result) {
    std::vector<ClusterDetailCsvRow> rows;
    rows.reserve(result.hierarchy.clusters.size());
    for (size_t cluster_index = 0; cluster_index < result.hierarchy.clusters.size(); cluster_index++) {
        const auto& cluster = result.hierarchy.clusters[cluster_index];
        ClusterDetailCsvRow row;
        row.cluster_id = cluster.id;
        row.level = cluster.level;
        row.active_detectors = cluster.active_detectors;
        row.diameter_bound = cluster.diameter_bound;
        row.buffer_bound = cluster.buffer_bound;
        row.diameter = cluster.diameter;
        row.nearest_boundary_match_distance = cluster.nearest_boundary_match_distance;
        row.nearest_external_detector_distance = cluster.nearest_external_detector_distance;
        row.forced_by_max_level = cluster.forced_by_max_level;
        row.event_count = value_or_zero(result.cluster_step_event_counts, cluster_index);
        row.critical_path_event_count = value_or_zero(result.cluster_critical_path_event_counts, cluster_index);
        row.has_cutoff_algorithmic_time = value_or_zero(result.cluster_has_cutoff_algorithmic_time, cluster_index);
        row.cutoff_algorithmic_time = value_or_zero(result.cluster_cutoff_algorithmic_times, cluster_index);
        row.prior_level_critical_event_count = value_or_zero(result.cluster_prior_level_critical_event_counts, cluster_index);
        row.overlap_event_base = value_or_zero(result.cluster_overlap_event_bases, cluster_index);
        row.events_at_or_before_cutoff = value_or_zero(result.cluster_events_at_or_before_cutoff, cluster_index);
        row.events_after_cutoff = value_or_zero(result.cluster_events_after_cutoff, cluster_index);
        row.stop_algorithmic_time = value_or_zero(result.cluster_stop_algorithmic_times, cluster_index);
        row.critical_path_algorithmic_time = value_or_zero(result.cluster_critical_path_algorithmic_times, cluster_index);
        rows.push_back(std::move(row));
    }
    return rows;
}

void append_mistake_log_rows_if_requested(
    const std::string& case_name,
    const CaseData& data,
    int rep_index,
    size_t shot_index,
    const std::string& method,
    uint64_t predicted_observable_mask,
    uint64_t expected_global_observable_mask,
    const pm::LockstepDecodeResult& result) {
    const char* shots_env = std::getenv("PYMATCHING_MISTAKE_SHOTS_CSV");
    const char* clusters_env = std::getenv("PYMATCHING_MISTAKE_CLUSTERS_CSV");
    const char* detectors_env = std::getenv("PYMATCHING_MISTAKE_DETECTORS_CSV");
    const char* matching_env = std::getenv("PYMATCHING_MISTAKE_MATCHING_CSV");
    const std::string shots_path = shots_env == nullptr ? "" : shots_env;
    const std::string clusters_path = clusters_env == nullptr ? "" : clusters_env;
    const std::string detectors_path = detectors_env == nullptr ? "" : detectors_env;
    const std::string matching_path = matching_env == nullptr ? "" : matching_env;
    if (shots_path.empty() && clusters_path.empty() && detectors_path.empty() && matching_path.empty()) {
        return;
    }
    const auto& shot = data.shots.at(shot_index);
    auto rows = make_cluster_detail_rows_from_result(result);
    if (!shots_path.empty()) {
        const bool write_header = !file_has_content(shots_path);
        std::ofstream out(shots_path, std::ios::app);
        if (!out) {
            throw std::runtime_error("failed to open mistake shots csv output: " + shots_path);
        }
        if (write_header) {
            out << "case,distance,rounds,noise,rep,shot,method,predicted_observable_mask,"
                   "expected_global_observable_mask,sample_obs_mask,active_detector_count,active_detectors,"
                   "cluster_count,max_level,parallel_cluster_events_total,parallel_level_critical_events,"
                   "parallel_nonempty_event_levels,parallel_critical_path_algorithmic_time,root_obs_crossed\n";
        }
        size_t max_level = 0;
        uint64_t total_cluster_events = 0;
        uint64_t max_critical_events = 0;
        for (const auto& row : rows) {
            max_level = std::max(max_level, row.level);
            total_cluster_events += row.event_count;
            max_critical_events = std::max(max_critical_events, row.critical_path_event_count);
        }
        uint64_t nonempty_levels = 0;
        std::set<size_t> levels_with_events;
        for (const auto& row : rows) {
            if (row.event_count > 0) levels_with_events.insert(row.level);
        }
        nonempty_levels = levels_with_events.size();
        std::vector<uint64_t> obs_bits;
        for (size_t k = 0; k < result.root_aggregate_result.obs_crossed.size(); k++) {
            if (result.root_aggregate_result.obs_crossed[k]) obs_bits.push_back(k);
        }
        out << case_name
            << "," << data.distance
            << "," << data.rounds
            << "," << data.noise
            << "," << rep_index
            << "," << shot_index
            << "," << method
            << "," << predicted_observable_mask
            << "," << expected_global_observable_mask
            << "," << shot.obs_mask
            << "," << shot.hits.size()
            << "," << join_u64_semicolon(shot.hits)
            << "," << rows.size()
            << "," << max_level
            << "," << total_cluster_events
            << "," << max_critical_events
            << "," << nonempty_levels
            << "," << value_or_zero(result.cluster_critical_path_algorithmic_times, rows.empty() ? 0 : rows.front().cluster_id)
            << "," << join_u64_semicolon(obs_bits)
            << "\n";
    }
    append_cluster_detail_csv_rows(clusters_path, detectors_path, "", method, case_name, data, rep_index, shot_index, rows);
    if (!matching_path.empty()) {
        const bool write_header = !file_has_content(matching_path);
        std::ofstream out(matching_path, std::ios::app);
        if (!out) {
            throw std::runtime_error("failed to open mistake matching csv output: " + matching_path);
        }
        if (write_header) {
            out << "case,distance,rounds,noise,rep,shot,method,cluster_id,level,execution_state,"
                   "active_detector_count,active_detectors,extraction_detector_count,extraction_detectors,"
                   "cluster_result_is_provisional,has_stopped_state,stopped_algorithmic_time,"
                   "stopped_detector_state_count,stopped_region_state_count,stopped_top_region_count,"
                   "stopped_region_match_count,stopped_blossom_child_edge_count,"
                   "stopped_negative_weight_obs_mask,stopped_negative_weight_obs_bits,"
                   "has_debug_compact_result,debug_compact_obs_mask,debug_compact_obs_bits,debug_compact_weight,"
                   "has_direct_compact_result,direct_compact_obs_mask,direct_compact_obs_bits,direct_compact_weight,"
                   "has_direct_result,direct_result_obs_bits,direct_result_weight,"
                   "stopped_region_matches,stopped_detector_ownership\n";
        }
        for (size_t cluster_index = 0; cluster_index < result.hierarchy.clusters.size(); cluster_index++) {
            const auto& cluster = result.hierarchy.clusters[cluster_index];
            const bool has_stopped = cluster_index < result.cluster_has_stopped_state.size() && result.cluster_has_stopped_state[cluster_index];
            const auto* state = has_stopped && cluster_index < result.cluster_stopped_states.size()
                                  ? &result.cluster_stopped_states[cluster_index]
                                  : nullptr;
            const bool has_debug = cluster_index < result.cluster_has_debug_compact_result.size() &&
                                   result.cluster_has_debug_compact_result[cluster_index];
            const bool has_direct_compact = cluster_index < result.cluster_has_direct_compact_result.size() &&
                                            result.cluster_has_direct_compact_result[cluster_index];
            const bool has_direct = cluster_index < result.cluster_has_direct_result.size() &&
                                    result.cluster_has_direct_result[cluster_index];
            const auto debug_result = has_debug ? result.cluster_debug_compact_results[cluster_index] : pm::MatchingResult();
            const auto direct_compact_result = has_direct_compact ? result.cluster_direct_compact_results[cluster_index] : pm::MatchingResult();
            out << case_name
                << "," << data.distance
                << "," << data.rounds
                << "," << data.noise
                << "," << rep_index
                << "," << shot_index
                << "," << method
                << "," << cluster.id
                << "," << cluster.level
                << "," << (cluster_index < result.cluster_execution_states.size()
                              ? lockstep_state_name(result.cluster_execution_states[cluster_index])
                              : std::string("UNKNOWN"))
                << "," << cluster.active_detectors.size()
                << "," << join_u64_semicolon(cluster.active_detectors)
                << "," << (cluster_index < result.cluster_extraction_detection_events.size()
                              ? result.cluster_extraction_detection_events[cluster_index].size()
                              : 0)
                << "," << (cluster_index < result.cluster_extraction_detection_events.size()
                              ? join_u64_semicolon(result.cluster_extraction_detection_events[cluster_index])
                              : std::string())
                << "," << (cluster_index < result.cluster_result_is_provisional.size() && result.cluster_result_is_provisional[cluster_index] ? 1 : 0)
                << "," << (has_stopped ? 1 : 0)
                << "," << (state == nullptr ? 0 : state->algorithmic_time)
                << "," << (state == nullptr ? 0 : state->detector_states.size())
                << "," << (state == nullptr ? 0 : state->region_states.size())
                << "," << (state == nullptr ? 0 : stopped_top_region_count(*state))
                << "," << (state == nullptr ? 0 : stopped_region_match_count(*state))
                << "," << (state == nullptr ? 0 : stopped_blossom_child_edge_count(*state))
                << "," << (state == nullptr ? 0 : (uint64_t)state->negative_weight_obs_mask)
                << "," << (state == nullptr ? std::string() : join_obs_bits_from_mask(state->negative_weight_obs_mask))
                << "," << (has_debug ? 1 : 0)
                << "," << (has_debug ? (uint64_t)debug_result.obs_mask : 0)
                << "," << (has_debug ? join_obs_bits_from_mask(debug_result.obs_mask) : std::string())
                << "," << (has_debug ? debug_result.weight : 0)
                << "," << (has_direct_compact ? 1 : 0)
                << "," << (has_direct_compact ? (uint64_t)direct_compact_result.obs_mask : 0)
                << "," << (has_direct_compact ? join_obs_bits_from_mask(direct_compact_result.obs_mask) : std::string())
                << "," << (has_direct_compact ? direct_compact_result.weight : 0)
                << "," << (has_direct ? 1 : 0)
                << "," << (has_direct ? join_bool_obs_bits(result.cluster_direct_results[cluster_index].obs_crossed) : std::string())
                << "," << (has_direct ? result.cluster_direct_results[cluster_index].weight : 0)
                << "," << (state == nullptr ? std::string() : stopped_region_match_summary(*state))
                << "," << (state == nullptr ? std::string() : stopped_detector_ownership_summary(*state))
                << "\n";
        }
    }
}

void append_fixed_cluster_detail_csv_rows_if_requested(
    const std::string& case_name,
    const CaseData& data,
    int rep_index,
    size_t shot_index,
    const pm::LockstepDecodeResult& result) {
    const char* members_env = std::getenv("PYMATCHING_CLUSTER_MEMBERS_CSV");
    const char* assignments_env = std::getenv("PYMATCHING_DETECTOR_ASSIGNMENTS_CSV");
    const char* level_params_env = std::getenv("PYMATCHING_LEVEL_PARAMS_CSV");
    const std::string members_path = members_env == nullptr ? "" : members_env;
    const std::string assignments_path = assignments_env == nullptr ? "" : assignments_env;
    const std::string level_params_path = level_params_env == nullptr ? "" : level_params_env;
    if (members_path.empty() && assignments_path.empty() && level_params_path.empty()) {
        return;
    }
    auto rows = make_cluster_detail_rows_from_result(result);
    append_cluster_detail_csv_rows(
        members_path, assignments_path, level_params_path,
        method_label_from_env("fixed"), case_name, data, rep_index, shot_index, rows);
}

uint64_t ns_since(Clock::time_point start) {
    if (eventcount_wall_timing_disabled()) {
        return 0;
    }
    return (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count();
}

CaseData read_case(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("failed to open " + path);
    size_t num_nodes, num_obs, num_edges, num_shots, rounds, distance;
    double noise;
    f >> num_nodes >> num_obs >> num_edges >> num_shots >> rounds >> distance >> noise;
    pm::UserGraph graph(num_nodes, num_obs);
    for (size_t k = 0; k < num_edges; k++) {
        char tag; int64_t u, v; double w, p; size_t nf;
        f >> tag >> u >> v >> w >> p >> nf;
        std::vector<size_t> obs(nf);
        for (size_t j = 0; j < nf; j++) f >> obs[j];
        if (tag != 'e') throw std::runtime_error("bad edge tag");
        if (v < 0) graph.add_or_merge_boundary_edge((size_t)u, obs, w, p, pm::INDEPENDENT);
        else graph.add_or_merge_edge((size_t)u, (size_t)v, obs, w, p, pm::INDEPENDENT);
    }
    graph.loaded_from_dem_without_correlations = true;
    std::vector<Shot> shots;
    shots.reserve(num_shots);
    for (size_t k = 0; k < num_shots; k++) {
        char tag; uint64_t obs; size_t nh;
        f >> tag >> obs >> nh;
        if (tag != 's') throw std::runtime_error("bad shot tag");
        Shot shot; shot.obs_mask = obs; shot.hits.resize(nh);
        for (size_t j = 0; j < nh; j++) f >> shot.hits[j];
        shots.push_back(std::move(shot));
    }
    return {std::move(graph), std::move(shots), rounds, distance, noise, num_edges};
}

uint64_t saturated_subtract(uint64_t a, uint64_t b) { return a > b ? a - b : 0; }
uint64_t max_level_cluster_work_ns(const pm::LockstepProfilingStats& stats) {
    return std::max(stats.max_level_cluster_step_wall_ns, stats.single_root_decode_wall_ns);
}
uint64_t max_level_cluster_lifecycle_ns(const pm::LockstepProfilingStats& stats) {
    return std::max(stats.max_level_cluster_lifecycle_wall_ns, stats.single_root_decode_wall_ns);
}

// Canonical comparison metric for future measurements.
//
// Execution model: every cluster at every level may run concurrently, and each
// worker already owns its detector graph/MWPM state before online sparse-blossom
// execution starts.  Therefore the critical path is only the sparse-blossom
// execution section on the slowest cluster over all levels, including child-state
// imports/checkpoints/export, but excluding detector-graph distribution,
// subgraph/materialization work, and initial MWPM setup.
uint64_t parallel_policy_critical_path_lifecycle_ns(const pm::LockstepProfilingStats& stats) {
    return max_level_cluster_lifecycle_ns(stats);
}

// Runtime-cost stage split for optimizing the decoder.
//
// The parallel decoder cost is split into:
//   1. precompute: graph/circuit data that can be prepared once and reused,
//      currently measured as clustering_precomputable_distance_lookup_wall_ns;
//   2. preprocessing: per-shot work before/around the sparse-blossom critical
//      path, excluding the reusable precompute;
//   3. parallel sparse blossom: the policy critical path, i.e. the slowest
//      maximum-level cluster execution time under all-level/all-cluster
//      concurrency. This includes imports/checkpoints/export but excludes
//      initial setup and pre-distributed worker graph state.
//
// For runtime optimization, use
// runtime_preprocessing_plus_parallel_sparse_blossom_ms_per_shot. It intentionally
// excludes the reusable precompute and is the sum of preprocessing plus the
// parallel sparse-blossom critical-path lifecycle.
uint64_t runtime_preprocessing_excluding_precompute_ns(const pm::LockstepProfilingStats& stats) {
    return saturated_subtract(stats.clustering_wall_ns, stats.clustering_precomputable_distance_lookup_wall_ns) +
        stats.root_extraction_wall_ns;
}
uint64_t runtime_preprocessing_plus_parallel_sparse_blossom_ns(const pm::LockstepProfilingStats& stats) {
    return runtime_preprocessing_excluding_precompute_ns(stats) + parallel_policy_critical_path_lifecycle_ns(stats);
}
// User-facing online ideal time for recurring per-shot work.
//
// This is the quantity to use when asking "how long from this shot's active
// detector pattern until the observable prediction is available?"  It excludes
// graph/circuit work that can be prepared once and reused, such as detector-
// graph distribution and graph-distance/radius-neighbor precomputation, but it
// does include all shot-dependent preprocessing plus the policy sparse-blossom
// critical path.
uint64_t ideal_parallel_time_ns(const pm::LockstepProfilingStats& stats) {
    return runtime_preprocessing_plus_parallel_sparse_blossom_ns(stats);
}
uint64_t ideal_runtime_preprocessing_plus_parallel_time_ns(const pm::LockstepProfilingStats& stats) {
    return runtime_preprocessing_plus_parallel_sparse_blossom_ns(stats);
}
uint64_t scheduler_bookkeeping_wall_ns(const pm::LockstepProfilingStats& stats) {
    uint64_t tracked = stats.cluster_step_wall_ns + stats.child_import_wall_ns +
        stats.overlap_checkpoint_capture_wall_ns + stats.overlap_checkpoint_restore_wall_ns;
    return saturated_subtract(stats.scheduler_wall_ns, tracked);
}
uint64_t parallel_overhead_ns(const pm::LockstepProfilingStats& stats) {
    return stats.clustering_wall_ns + stats.subgraph_build_wall_ns + stats.initial_mwpm_setup_wall_ns +
        stats.root_extraction_wall_ns + stats.child_import_wall_ns + stats.overlap_checkpoint_capture_wall_ns +
        stats.overlap_checkpoint_restore_wall_ns + scheduler_bookkeeping_wall_ns(stats);
}
uint64_t ideal_parallel_ns(const pm::LockstepProfilingStats& stats) {
    return parallel_overhead_ns(stats) + max_level_cluster_work_ns(stats);
}
uint64_t ideal_parallel_excluding_precompute_ns(const pm::LockstepProfilingStats& stats) {
    return saturated_subtract(ideal_parallel_ns(stats), stats.clustering_precomputable_distance_lookup_wall_ns);
}
uint64_t overhead_excluding_precompute_ns(const pm::LockstepProfilingStats& stats) {
    return saturated_subtract(parallel_overhead_ns(stats), stats.clustering_precomputable_distance_lookup_wall_ns);
}
uint64_t single_core_preprocessing_excluding_precompute_ns(const pm::LockstepProfilingStats& stats) {
    uint64_t preprocessing = stats.clustering_wall_ns + stats.subgraph_build_wall_ns +
        stats.initial_mwpm_setup_wall_ns + stats.root_extraction_wall_ns +
        stats.overlap_checkpoint_capture_wall_ns + stats.overlap_checkpoint_restore_wall_ns;
    return saturated_subtract(preprocessing, stats.clustering_precomputable_distance_lookup_wall_ns);
}
uint64_t multi_worker_ideal_excluding_precompute_with_comm_ns(const pm::LockstepProfilingStats& stats) {
    return single_core_preprocessing_excluding_precompute_ns(stats) +
        stats.ideal_worker_cluster_makespan_wall_ns +
        stats.ideal_worker_communication_estimate_ns;
}

struct Totals {
    uint64_t wall_ns = 0;
    uint64_t parallel_policy_critical_path_lifecycle_ns = 0;
    uint64_t ideal_parallel_time_ns = 0;
    uint64_t ideal_runtime_preprocessing_plus_parallel_time_ns = 0;
    uint64_t runtime_preprocessing_plus_parallel_sparse_blossom_ns = 0;
    uint64_t runtime_preprocessing_excluding_precompute_ns = 0;
    uint64_t runtime_parallel_sparse_blossom_ns = 0;
    uint64_t clustering_ns = 0;
    uint64_t component_ns = 0;
    uint64_t diameter_ns = 0;
    uint64_t precomputable_distance_lookup_ns = 0;
    uint64_t influence_ns = 0;
    uint64_t parent_ns = 0;
    uint64_t used_active_pair_table = 0;
    uint64_t used_sparse_frontier = 0;
    uint64_t used_component_lookup = 0;
    uint64_t forced_max_level_shortcut_count = 0;
    uint64_t processed_frontier_events = 0;
    uint64_t generated_collision_events = 0;
    uint64_t processed_collision_events = 0;
    uint64_t exact_internal_pair_checks = 0;
    uint64_t exact_external_pair_checks = 0;
    uint64_t exact_diameter_checks = 0;
    uint64_t max_component_active_size = 0;
    uint64_t max_residual_active_size = 0;
    uint64_t subgraph_build_ns = 0;
    uint64_t worker_graph_distribution_excluded_ns = 0;
    uint64_t initial_setup_ns = 0;
    uint64_t scheduler_ns = 0;
    uint64_t child_import_ns = 0;
    uint64_t checkpoint_capture_ns = 0;
    uint64_t checkpoint_restore_ns = 0;
    uint64_t cluster_step_ns = 0;
    uint64_t max_level_cluster_ns = 0;
    uint64_t max_level_lifecycle_ns = 0;
    uint64_t max_level_lifecycle_initial_setup_ns = 0;
    uint64_t max_level_lifecycle_child_import_ns = 0;
    uint64_t max_level_lifecycle_checkpoint_capture_ns = 0;
    uint64_t max_level_lifecycle_checkpoint_restore_ns = 0;
    uint64_t max_level_lifecycle_mark_stop_ns = 0;
    uint64_t max_level_lifecycle_cluster_id = 0;
    uint64_t overhead_ns = 0;
    uint64_t ideal_ns = 0;
    uint64_t overhead_excluding_precompute_ns = 0;
    uint64_t ideal_excluding_precompute_ns = 0;
    uint64_t single_core_preprocessing_excluding_precompute_ns = 0;
    uint64_t ideal_worker_cluster_makespan_ns = 0;
    uint64_t ideal_worker_communication_estimate_ns = 0;
    uint64_t multi_worker_ideal_excluding_precompute_with_comm_ns = 0;
    uint64_t ideal_worker_count_used = 0;
    uint64_t child_import_count = 0;
    uint64_t clusters = 0;
    uint64_t max_level = 0;
    uint64_t peak_width = 0;
    uint64_t max_cluster_size = 0;
    uint64_t scheduler_batches = 0;
    uint64_t scheduler_cluster_steps = 0;
    uint64_t direct_worker_init_ns = 0;
    uint64_t direct_worker_advance_ns = 0;
    uint64_t direct_worker_extract_ns = 0;
    uint64_t direct_worker_release_ns = 0;
    uint64_t direct_worker_total_ns = 0;
    uint64_t max_direct_worker_init_ns = 0;
    uint64_t max_direct_worker_advance_ns = 0;
    uint64_t max_direct_worker_extract_ns = 0;
    uint64_t max_direct_worker_release_ns = 0;
    uint64_t max_direct_worker_total_ns = 0;
    uint64_t direct_worker_active_detectors = 0;
    uint64_t max_direct_worker_active_detectors = 0;
    uint64_t global_mwpm_events = 0;
    long double global_mwpm_events_square_sum = 0;
    uint64_t parallel_cluster_events = 0;
    long double parallel_cluster_events_square_sum = 0;
    uint64_t parallel_level_max_event_sum = 0;
    long double parallel_level_max_event_square_sum = 0;
    uint64_t parallel_max_cluster_events = 0;
    uint64_t parallel_nonempty_event_levels = 0;
    uint64_t parallel_critical_path_algorithmic_time = 0;
    long double parallel_critical_path_algorithmic_time_square_sum = 0;
    size_t mistakes = 0;
    size_t exceptions = 0;
    uint64_t warmup_wall_ns = 0;
    uint64_t execution_wall_ns = 0;
};


void measure_global_stage_breakdown(pm::UserGraph& graph, const std::vector<Shot>& shots, int reps, Totals& t) {
    (void)graph.get_mwpm();
    for (int r = 0; r < reps; r++) {
        for (const auto& shot : shots) {
            auto& mwpm = graph.get_mwpm();
            auto total_start = Clock::now();

            auto init_start = Clock::now();
            pm::initialize_mwpm_for_detection_events(mwpm, shot.hits);
            auto init_ns = ns_since(init_start);

            auto advance_start = Clock::now();
            uint64_t event_count = 0;
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
                throw std::invalid_argument("No perfect matching in global stage breakdown.");
            }
            auto advance_ns = ns_since(advance_start);

            auto extract_start = Clock::now();
            (void)pm::extract_compact_result_from_current_mwpm_state(mwpm, shot.hits);
            auto extract_ns = ns_since(extract_start);

            auto total_ns = ns_since(total_start);
            t.direct_worker_init_ns += init_ns;
            t.direct_worker_advance_ns += advance_ns;
            t.direct_worker_extract_ns += extract_ns;
            t.direct_worker_total_ns += total_ns;
            t.max_direct_worker_init_ns += init_ns;
            t.max_direct_worker_advance_ns += advance_ns;
            t.max_direct_worker_extract_ns += extract_ns;
            t.max_direct_worker_total_ns += total_ns;
            t.direct_worker_active_detectors += shot.hits.size();
            t.max_direct_worker_active_detectors += shot.hits.size();
            t.global_mwpm_events += event_count;
            t.global_mwpm_events_square_sum += static_cast<long double>(event_count) * static_cast<long double>(event_count);
        }
    }
}


void measure_global_event_counts_only(pm::UserGraph& graph, const std::vector<Shot>& shots, int reps, Totals& t) {
    (void)graph.get_mwpm();
    for (int r = 0; r < reps; r++) {
        for (const auto& shot : shots) {
            auto& mwpm = graph.get_mwpm();
            pm::initialize_mwpm_for_detection_events(mwpm, shot.hits);
            uint64_t event_count = 0;
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
                throw std::invalid_argument("No perfect matching in global event count.");
            }
            (void)pm::extract_compact_result_from_current_mwpm_state(mwpm, shot.hits);
            t.global_mwpm_events += event_count;
            t.global_mwpm_events_square_sum += static_cast<long double>(event_count) * static_cast<long double>(event_count);
        }
    }
}

pm::AdaptiveEventTiming count_mwpm_event_timing(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    pm::cumulative_time_int cutoff_algorithmic_time) {
    pm::initialize_mwpm_for_detection_events(mwpm, detection_events);
    pm::AdaptiveEventTiming timing;
    while (true) {
        auto event = mwpm.flooder.run_until_next_mwpm_notification();
        if (event.event_type == pm::NO_EVENT) {
            break;
        }
        timing.raw_event_count++;
        if (mwpm.flooder.queue.cur_time > cutoff_algorithmic_time) {
            timing.events_after_cutoff++;
        }
        mwpm.process_event(event);
    }
    if (mwpm.node_arena.allocated.size() != mwpm.node_arena.available.size()) {
        mwpm.reset();
        throw std::invalid_argument("No perfect matching in event-count-only MWPM run.");
    }
    timing.stop_algorithmic_time = mwpm.flooder.queue.cur_time;
    // Unlike the normal decode path, event-count-only does not extract a result.
    // Reset explicitly so the single reusable full-graph MWPM can be used for the next cluster.
    mwpm.reset();
    return timing;
}

uint64_t count_mwpm_events_only(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    return count_mwpm_event_timing(
        mwpm,
        detection_events,
        std::numeric_limits<pm::cumulative_time_int>::min()).raw_event_count;
}


Totals run_parallel_cluster_only(
    pm::UserGraph& graph,
    const std::vector<Shot>& shots,
    int reps,
    const pm::ProcessingClusterConfig& cluster_config) {
    Totals totals;
    const auto& matching_graph = graph.get_matching_graph_for_parallel_clustering();
    auto* cache = graph.get_processing_cluster_graph_cache();
    pm::ensure_processing_cluster_radius_neighbors_precomputed(matching_graph, *cache, cluster_config);
    const size_t max_level = pm::processing_cluster_effective_max_level(cluster_config);

    auto total_start = Clock::now();
    for (int r = 0; r < reps; r++) {
        for (const auto& shot : shots) {
            std::vector<uint64_t> residual = shot.hits;
            std::sort(residual.begin(), residual.end());
            residual.erase(std::unique(residual.begin(), residual.end()), residual.end());
            pm::ProcessingClusterConfig shot_cluster_config = cluster_config;

            std::vector<uint64_t> cluster_count_by_level(max_level + 2, 0);
            size_t shot_clusters = 0;
            size_t shot_max_level = 0;
            size_t shot_peak_width = 0;
            size_t shot_max_cluster_size = 0;
            size_t highest_non_root_level = 0;
            uint64_t shot_clustering_ns = 0;

            for (size_t level = 1; level <= max_level && !residual.empty(); level++) {
                auto start = eventcount_optional_now();
                auto result = pm::run_growing_only_clustering_with_adaptive_bounds(
                    matching_graph, std::move(residual), level, shot_cluster_config, cache);
                shot_clustering_ns += ns_since(start);
                for (const auto& cluster : result.accepted_clusters) {
                    if (level >= cluster_count_by_level.size()) {
                        cluster_count_by_level.resize(level + 1, 0);
                    }
                    cluster_count_by_level[level]++;
                    shot_clusters++;
                    highest_non_root_level = std::max(highest_non_root_level, level);
                    shot_max_level = std::max(shot_max_level, level);
                    shot_max_cluster_size = std::max(shot_max_cluster_size, cluster.active_detectors.size());
                }
                residual = std::move(result.residual_active_detectors);
            }

            const char* keep_empty_root_env = std::getenv("PYMATCHING_GROWING_ONLY_KEEP_EMPTY_TERMINAL_ROOT");
            const bool keep_empty_root = keep_empty_root_env != nullptr && std::string(keep_empty_root_env) == "1";
            if (shot_clusters == 0) {
                cluster_count_by_level[0]++;
                shot_clusters = 1;
                shot_max_level = 0;
                shot_max_cluster_size = residual.size();
            } else if (residual.empty() && !keep_empty_root) {
            } else {
                size_t root_level = std::max<size_t>(highest_non_root_level + 1, 1);
                if (root_level >= cluster_count_by_level.size()) {
                    cluster_count_by_level.resize(root_level + 1, 0);
                }
                cluster_count_by_level[root_level]++;
                shot_clusters++;
                shot_max_level = std::max(shot_max_level, root_level);
                shot_max_cluster_size = std::max(shot_max_cluster_size, residual.size());
            }

            for (size_t level = 0; level < cluster_count_by_level.size(); level++) {
                shot_peak_width = std::max<size_t>(shot_peak_width, cluster_count_by_level[level]);
            }
            totals.clustering_ns += shot_clustering_ns;
            totals.clusters += shot_clusters;
            totals.max_level += shot_max_level;
            totals.peak_width += shot_peak_width;
            totals.max_cluster_size += shot_max_cluster_size;
        }
    }
    totals.wall_ns = ns_since(total_start);
    totals.execution_wall_ns = totals.wall_ns;
    return totals;
}

Totals run_parallel_event_count_only(
    pm::UserGraph& graph,
    const std::vector<Shot>& shots,
    int reps,
    const pm::ProcessingClusterConfig& cluster_config) {
    Totals totals;
    const auto& matching_graph = graph.get_matching_graph_for_parallel_clustering();
    auto* cache = graph.get_processing_cluster_graph_cache();
    auto& mwpm = graph.get_mwpm();
    pm::ProcessingClusterConfig event_count_cluster_config = cluster_config;
    std::map<std::vector<uint64_t>, uint64_t> adaptive_event_count_cache;
    std::map<std::pair<std::vector<uint64_t>, pm::cumulative_time_int>, pm::AdaptiveEventTiming> adaptive_event_timing_cache;
    event_count_cluster_config.adaptive_event_count_callback = [&](const std::vector<uint64_t>& detectors) -> uint64_t {
        std::vector<uint64_t> key = detectors;
        std::sort(key.begin(), key.end());
        key.erase(std::unique(key.begin(), key.end()), key.end());
        auto it = adaptive_event_count_cache.find(key);
        if (it != adaptive_event_count_cache.end()) {
            return it->second;
        }
        uint64_t events = count_mwpm_events_only(mwpm, key);
        adaptive_event_count_cache.emplace(std::move(key), events);
        return events;
    };
    event_count_cluster_config.adaptive_event_timing_callback = [&](
        const std::vector<uint64_t>& detectors,
        pm::cumulative_time_int cutoff_algorithmic_time) -> pm::AdaptiveEventTiming {
        std::vector<uint64_t> key = detectors;
        std::sort(key.begin(), key.end());
        key.erase(std::unique(key.begin(), key.end()), key.end());
        auto cache_key = std::make_pair(key, cutoff_algorithmic_time);
        auto it = adaptive_event_timing_cache.find(cache_key);
        if (it != adaptive_event_timing_cache.end()) {
            return it->second;
        }
        pm::AdaptiveEventTiming timing = count_mwpm_event_timing(mwpm, key, cutoff_algorithmic_time);
        adaptive_event_timing_cache.emplace(std::move(cache_key), timing);
        adaptive_event_count_cache.emplace(key, timing.raw_event_count);
        return timing;
    };
    pm::ensure_processing_cluster_radius_neighbors_precomputed(matching_graph, *cache, event_count_cluster_config);
    const size_t max_level = pm::processing_cluster_effective_max_level(event_count_cluster_config);

    for (int r = 0; r < reps; r++) {
        for (const auto& shot : shots) {
            std::vector<uint64_t> residual = shot.hits;
            std::sort(residual.begin(), residual.end());
            residual.erase(std::unique(residual.begin(), residual.end()), residual.end());
            pm::ProcessingClusterConfig shot_cluster_config = event_count_cluster_config;

            std::vector<uint64_t> max_events_by_level(max_level + 2, 0);
            std::vector<uint64_t> cluster_count_by_level(max_level + 2, 0);
            size_t shot_clusters = 0;
            size_t shot_max_level = 0;
            size_t shot_peak_width = 0;
            size_t shot_max_cluster_size = 0;
            size_t highest_non_root_level = 0;

            for (size_t level = 1; level <= max_level && !residual.empty(); level++) {
                auto result = pm::run_growing_only_clustering_with_adaptive_bounds(
                    matching_graph, std::move(residual), level, shot_cluster_config, cache);
                for (const auto& cluster : result.accepted_clusters) {
                    const uint64_t event_count = count_mwpm_events_only(mwpm, cluster.active_detectors);
                    totals.parallel_cluster_events += event_count;
                    totals.parallel_max_cluster_events = std::max(totals.parallel_max_cluster_events, event_count);
                    if (level >= max_events_by_level.size()) {
                        max_events_by_level.resize(level + 1, 0);
                        cluster_count_by_level.resize(level + 1, 0);
                    }
                    max_events_by_level[level] = std::max(max_events_by_level[level], event_count);
                    cluster_count_by_level[level]++;
                    shot_clusters++;
                    highest_non_root_level = std::max(highest_non_root_level, level);
                    shot_max_level = std::max(shot_max_level, level);
                    shot_max_cluster_size = std::max(shot_max_cluster_size, cluster.active_detectors.size());
                }
                residual = std::move(result.residual_active_detectors);
            }

            const char* keep_empty_root_env = std::getenv("PYMATCHING_GROWING_ONLY_KEEP_EMPTY_TERMINAL_ROOT");
            const bool keep_empty_root = keep_empty_root_env != nullptr && std::string(keep_empty_root_env) == "1";
            if (shot_clusters == 0) {
                const uint64_t event_count = count_mwpm_events_only(mwpm, residual);
                totals.parallel_cluster_events += event_count;
                totals.parallel_max_cluster_events = std::max(totals.parallel_max_cluster_events, event_count);
                max_events_by_level[0] = std::max(max_events_by_level[0], event_count);
                cluster_count_by_level[0]++;
                shot_clusters = 1;
                shot_max_level = 0;
                shot_max_cluster_size = residual.size();
            } else if (residual.empty() && !keep_empty_root) {
            } else {
                size_t root_level = std::max<size_t>(highest_non_root_level + 1, 1);
                if (root_level >= max_events_by_level.size()) {
                    max_events_by_level.resize(root_level + 1, 0);
                    cluster_count_by_level.resize(root_level + 1, 0);
                }
                const uint64_t event_count = count_mwpm_events_only(mwpm, residual);
                totals.parallel_cluster_events += event_count;
                totals.parallel_max_cluster_events = std::max(totals.parallel_max_cluster_events, event_count);
                max_events_by_level[root_level] = std::max(max_events_by_level[root_level], event_count);
                cluster_count_by_level[root_level]++;
                shot_clusters++;
                shot_max_level = std::max(shot_max_level, root_level);
                shot_max_cluster_size = std::max(shot_max_cluster_size, residual.size());
            }

            uint64_t shot_max_cluster_events = 0;
            uint64_t shot_nonempty_event_levels = 0;
            for (size_t level = 0; level < max_events_by_level.size(); level++) {
                shot_max_cluster_events = std::max(shot_max_cluster_events, max_events_by_level[level]);
                if (max_events_by_level[level] > 0) {
                    shot_nonempty_event_levels++;
                }
                shot_peak_width = std::max<size_t>(shot_peak_width, cluster_count_by_level[level]);
            }
            totals.parallel_level_max_event_sum += shot_max_cluster_events;
            totals.parallel_nonempty_event_levels += shot_nonempty_event_levels;
            totals.clusters += shot_clusters;
            totals.max_level += shot_max_level;
            totals.peak_width += shot_peak_width;
            totals.max_cluster_size += shot_max_cluster_size;
        }
    }
    return totals;
}

uint64_t run_global(pm::UserGraph& graph, const std::vector<Shot>& shots, int reps, size_t& mistakes) {
    (void)graph.get_mwpm();
    uint64_t total = 0;
    std::vector<uint8_t> obs(graph.get_num_observables());
    for (int r = 0; r < reps; r++) {
        for (const auto& shot : shots) {
            std::fill(obs.begin(), obs.end(), 0);
            pm::total_weight_int weight = 0;
            auto start = eventcount_optional_now();
            pm::decode_detection_events(graph.get_mwpm(), shot.hits, obs.data(), weight, false);
            total += ns_since(start);
            uint64_t pred = 0;
            for (size_t k = 0; k < obs.size() && k < 64; k++) pred ^= (uint64_t)obs[k] << k;
            if (pred != shot.obs_mask) mistakes++;
        }
    }
    return total;
}

// Computes only the predicted observable flip bitmask from the global decoder.
// This is intentionally not a comparison of the full MWPM edge/matching configuration.
std::vector<uint64_t> compute_global_observable_predictions(pm::UserGraph& graph, const std::vector<Shot>& shots) {
    (void)graph.get_mwpm();
    std::vector<uint64_t> predictions;
    predictions.reserve(shots.size());
    std::vector<uint8_t> obs(graph.get_num_observables());
    for (const auto& shot : shots) {
        std::fill(obs.begin(), obs.end(), 0);
        pm::total_weight_int weight = 0;
        pm::decode_detection_events(graph.get_mwpm(), shot.hits, obs.data(), weight, false);
        uint64_t pred = 0;
        for (size_t k = 0; k < obs.size() && k < 64; k++) pred ^= (uint64_t)obs[k] << k;
        predictions.push_back(pred);
    }
    return predictions;
}

std::vector<std::string> split_nonempty(const std::string& text, char delimiter) {
    std::vector<std::string> parts;
    std::stringstream ss(text);
    std::string part;
    while (std::getline(ss, part, delimiter)) {
        if (!part.empty()) {
            parts.push_back(part);
        }
    }
    return parts;
}

std::vector<pm::ProcessingClusterBounds> parse_explicit_cluster_bounds(const std::string& spec) {
    std::vector<pm::ProcessingClusterBounds> bounds;
    for (const auto& item : split_nonempty(spec, ',')) {
        auto colon = item.find(':');
        if (colon == std::string::npos) {
            throw std::runtime_error(
                "PYMATCHING_CLUSTER_BOUNDS entries must have the form d:b, e.g. 10:25,20:60.");
        }
        auto d = static_cast<pm::cumulative_time_int>(std::stoll(item.substr(0, colon)));
        auto b = static_cast<pm::cumulative_time_int>(std::stoll(item.substr(colon + 1)));
        bounds.push_back(pm::ProcessingClusterBounds{d, b});
    }
    pm::ProcessingClusterConfig check_config;
    check_config.max_level = bounds.size();
    check_config.explicit_bounds_by_level = bounds;
    pm::validate_processing_cluster_stopping_lemma_bounds(check_config);
    return bounds;
}

pm::cumulative_time_int ceil_positive_scaled_distance(
    pm::cumulative_time_int unit,
    long double ratio) {
    if (!(unit > 0)) {
        throw std::runtime_error("linear diameter candidate unit must be positive");
    }
    if (!(ratio > 0) || !std::isfinite(static_cast<double>(ratio))) {
        throw std::runtime_error("linear diameter candidate ratio must be positive and finite");
    }
    long double scaled = static_cast<long double>(unit) * ratio;
    auto max_distance = static_cast<long double>(std::numeric_limits<pm::cumulative_time_int>::max() / 8);
    if (!(scaled > 0) || scaled >= max_distance) {
        throw std::runtime_error("linear diameter candidate start would overflow cumulative_time_int");
    }
    long double rounded = std::round(scaled);
    long double tolerance = std::max(1e-9L, std::fabs(scaled) * 1e-15L);
    if (std::fabs(scaled - rounded) <= tolerance) {
        scaled = rounded;
    }
    return static_cast<pm::cumulative_time_int>(std::ceil(scaled));
}

pm::cumulative_time_int detector_graph_max_edge_weight_plus_one(pm::UserGraph& graph) {
    const auto& matching_graph = graph.get_matching_graph_for_parallel_clustering();
    pm::cumulative_time_int max_weight = 0;
    for (const auto& node : matching_graph.nodes) {
        for (auto weight : node.neighbor_weights) {
            max_weight = std::max(max_weight, static_cast<pm::cumulative_time_int>(weight));
        }
    }
    if (!(max_weight > 0)) {
        throw std::runtime_error(
            "Could not set d1=w_max+1 because the detector graph has no positive edge weight.");
    }
    const auto max_distance = std::numeric_limits<pm::cumulative_time_int>::max() / 8;
    if (max_weight >= max_distance - 1) {
        throw std::runtime_error("d1=w_max+1 would overflow cumulative_time_int.");
    }
    return max_weight + 1;
}

std::vector<pm::cumulative_time_int> make_linear_diameter_candidates(
    pm::cumulative_time_int origin_unit,
    long double start_ratio,
    size_t candidate_count) {
    if (!(origin_unit > 0)) {
        throw std::runtime_error("linear diameter candidate origin unit must be positive");
    }
    auto start = ceil_positive_scaled_distance(origin_unit, start_ratio);
    std::vector<pm::cumulative_time_int> candidates;
    candidates.reserve(candidate_count);
    auto max_distance = std::numeric_limits<pm::cumulative_time_int>::max() / 8;
    for (size_t k = 0; k < candidate_count; k++) {
        if (k > static_cast<size_t>((max_distance - start) / origin_unit)) {
            break;
        }
        candidates.push_back(static_cast<pm::cumulative_time_int>(start + k * static_cast<size_t>(origin_unit)));
    }
    return candidates;
}

struct MatchingPathDistanceStats {
    size_t shots_seen = 0;
    size_t shots_with_hits = 0;
    size_t total_match_edges = 0;
    size_t detector_pair_match_edges = 0;
    size_t boundary_match_edges = 0;
    std::vector<uint64_t> detector_pair_distances;
    std::vector<uint64_t> all_match_distances;
};

double mean_uint64(const std::vector<uint64_t>& values) {
    if (values.empty()) {
        return 0;
    }
    long double sum = 0;
    for (auto v : values) {
        sum += static_cast<long double>(v);
    }
    return static_cast<double>(sum / static_cast<long double>(values.size()));
}

uint64_t percentile_uint64(std::vector<uint64_t> values, double p) {
    if (values.empty()) {
        return 0;
    }
    std::sort(values.begin(), values.end());
    size_t index = static_cast<size_t>(std::floor(p * static_cast<double>(values.size() - 1)));
    return values[index];
}

void print_matching_path_stats(const std::string& case_name, const MatchingPathDistanceStats& stats) {
    auto print_one = [&](const char* label, const std::vector<uint64_t>& values) {
        std::cerr << "MATCH_PATH_STATS case=" << case_name
                  << " kind=" << label
                  << " shots=" << stats.shots_seen
                  << " shots_with_hits=" << stats.shots_with_hits
                  << " count=" << values.size();
        if (values.empty()) {
            std::cerr << " mean=nan min=nan p50=nan p90=nan p99=nan max=nan\n";
            return;
        }
        auto sorted = values;
        std::sort(sorted.begin(), sorted.end());
        std::cerr << " mean=" << std::setprecision(17) << mean_uint64(values)
                  << " min=" << sorted.front()
                  << " p50=" << percentile_uint64(sorted, 0.50)
                  << " p90=" << percentile_uint64(sorted, 0.90)
                  << " p99=" << percentile_uint64(sorted, 0.99)
                  << " max=" << sorted.back()
                  << "\n";
    };
    std::cerr << "MATCH_PATH_COUNTS case=" << case_name
              << " shots=" << stats.shots_seen
              << " total_match_edges=" << stats.total_match_edges
              << " detector_pair_match_edges=" << stats.detector_pair_match_edges
              << " boundary_match_edges=" << stats.boundary_match_edges << "\n";
    print_one("detector_pair", stats.detector_pair_distances);
    print_one("all", stats.all_match_distances);
}

MatchingPathDistanceStats measure_matching_path_distances(
    pm::UserGraph& graph,
    const std::vector<Shot>& shots,
    size_t max_shots) {
    MatchingPathDistanceStats stats;
    auto& mwpm = graph.get_mwpm_with_search_graph();
    size_t limit = max_shots == 0 ? shots.size() : std::min(max_shots, shots.size());
    for (size_t shot_index = 0; shot_index < limit; shot_index++) {
        const auto& shot = shots[shot_index];
        stats.shots_seen++;
        if (!shot.hits.empty()) {
            stats.shots_with_hits++;
        }
        pm::decode_detection_events_to_match_edges(mwpm, shot.hits);
        stats.total_match_edges += mwpm.flooder.match_edges.size();
        for (const auto& match_edge : mwpm.flooder.match_edges) {
            size_t from = static_cast<size_t>(match_edge.loc_from - &mwpm.flooder.graph.nodes[0]);
            size_t to = match_edge.loc_to ? static_cast<size_t>(match_edge.loc_to - &mwpm.flooder.graph.nodes[0]) : SIZE_MAX;
            uint64_t path_weight = 0;
            mwpm.search_flooder.iter_edges_on_shortest_path_from_middle(
                from, to, [&](const pm::SearchGraphEdge& edge) {
                    path_weight += static_cast<uint64_t>(edge.detector_node->neighbor_weights[edge.neighbor_index]);
                });
            stats.all_match_distances.push_back(path_weight);
            if (match_edge.loc_to) {
                stats.detector_pair_match_edges++;
                stats.detector_pair_distances.push_back(path_weight);
            } else {
                stats.boundary_match_edges++;
            }
        }
    }
    return stats;
}

pm::cumulative_time_int choose_path_distance_unit(
    const MatchingPathDistanceStats& stats,
    const char* override_env) {
    if (override_env != nullptr) {
        return static_cast<pm::cumulative_time_int>(std::stoll(override_env));
    }
    const auto& preferred = !stats.detector_pair_distances.empty()
        ? stats.detector_pair_distances
        : stats.all_match_distances;
    double mean = mean_uint64(preferred);
    if (!(mean > 0)) {
        throw std::runtime_error(
            "Could not infer a path-distance unit because no matching paths were measured. "
            "Set PYMATCHING_CLUSTER_PATH_UNIT explicitly.");
    }
    return static_cast<pm::cumulative_time_int>(std::ceil(mean));
}

void print_selected_bounds(
    const std::string& case_name,
    pm::cumulative_time_int origin_unit,
    long double start_ratio,
    size_t candidate_count,
    double min_buffer_ratio,
    bool min_buffer_ratio_first_level_only,
    bool direct_min_d,
    bool optimize_final_d,
    bool phi_sequence,
    bool parameter_schedule,
    bool phi_floor_budget,
    double phi_floor,
    double phi_budget_per_level,
    const std::vector<pm::ProcessingClusterBounds>& bounds) {
    std::cerr << "PATH_SELECTED_BOUNDS case=" << case_name
              << " l_origin=" << origin_unit
              << " l_ratio=" << static_cast<double>(start_ratio)
              << " origin_unit=" << origin_unit
              << " start_ratio=" << static_cast<double>(start_ratio)
              << " start_d=" << (bounds.empty() ? 0 : bounds.front().diameter_bound)
              << " candidate_count=" << candidate_count
              << " selected_levels=" << bounds.size()
              << " min_buffer_ratio=" << min_buffer_ratio
              << " min_buffer_ratio_first_level_only="
              << (min_buffer_ratio_first_level_only ? 1 : 0)
              << " direct_min_d=" << (direct_min_d ? 1 : 0)
              << " optimize_final_d=" << (optimize_final_d ? 1 : 0)
              << " phi_sequence=" << (phi_sequence ? 1 : 0)
              << " parameter_schedule=" << (parameter_schedule ? 1 : 0)
              << " phi_floor_budget=" << (phi_floor_budget ? 1 : 0)
              << " phi_floor=" << phi_floor
              << " phi_budget_per_level=" << phi_budget_per_level << "\n";
    for (size_t k = 0; k < bounds.size(); k++) {
        const auto& b = bounds[k];
        std::cerr << "PATH_SELECTED_BOUND level=" << (k + 1)
                  << " d=" << b.diameter_bound
                  << " b=" << b.buffer_bound;
        if (b.has_long_double_schedule_bounds) {
            std::cerr << " schedule_d=" << static_cast<double>(b.schedule_diameter_bound)
                      << " schedule_b=" << static_cast<double>(b.schedule_buffer_bound)
                      << " runtime_clamped=" << (b.runtime_bounds_clamped ? 1 : 0);
        }
        if (b.has_phi_buffer_certificate) {
            std::cerr << " phi=" << b.phi
                      << " phi_required_b=" << b.phi_required_buffer_bound
                      << " phi_budget_required_b=" << b.phi_budget_required_buffer_bound
                      << " ratio_required_b=" << b.ratio_required_buffer_bound
                      << " min_required_b=" << b.min_required_buffer_bound
                      << " phi_cert_ok="
                      << (pm::processing_cluster_bounds_phi_buffer_certificate_ok(b) ? 1 : 0);
        } else {
            std::cerr << " phi_cert_ok=0";
        }
        if (k > 0) {
            std::cerr << " delta_d=" << (b.diameter_bound - bounds[k - 1].diameter_bound);
        }
        std::cerr << "\n";
    }
}


struct GrowingOnlyTotals {
    uint64_t event_generation_ns = 0;
    uint64_t event_processing_ns = 0;
    uint64_t wall_ns = 0;
    uint64_t generated_events = 0;
    uint64_t processed_events = 0;
    uint64_t internal_collisions = 0;
    uint64_t external_collisions = 0;
    uint64_t union_count = 0;
    uint64_t accepted_clusters = 0;
    uint64_t residual_active = 0;
    uint64_t removed_or_hidden_clusters = 0;
    uint64_t parity_rejected_clusters = 0;
    uint64_t diameter_rejected_clusters = 0;
    uint64_t max_active_rejected_clusters = 0;
    uint64_t forced_by_max_level_clusters = 0;
    uint64_t sparse_frontier_events = 0;
    uint64_t processed_frontier_events = 0;
    uint64_t generated_collision_events = 0;
    uint64_t processed_collision_events = 0;
    uint64_t exact_internal_pair_checks = 0;
    uint64_t exact_external_pair_checks = 0;
    uint64_t exact_diameter_checks = 0;
    uint64_t sparse_frontier_used = 0;
    uint64_t component_lookup_used = 0;
    uint64_t active_pair_table_used = 0;
    uint64_t input_active = 0;
    uint64_t accepted_active = 0;
    uint64_t accepted_cluster_size_sum = 0;
    uint64_t accepted_cluster_max_size = 0;
    pm::cumulative_time_int max_diameter = 0;
    pm::cumulative_time_int max_buffer = 0;
    pm::cumulative_time_int max_ready_time = 0;
    pm::cumulative_time_int max_first_external_collision_time = 0;
    size_t examples_printed = 0;
};

bool env_flag_enabled(const char* name) {
    const char* value = std::getenv(name);
    return value != nullptr && std::string(value) == "1";
}

void apply_adaptive_cluster_bounds_env(pm::ProcessingClusterConfig& config) {
    if (!env_flag_enabled("PYMATCHING_ADAPTIVE_CLUSTER_BOUNDS")) {
        return;
    }
    config.adaptive_cluster_bounds = true;
    if (const char* e = std::getenv("PYMATCHING_ADAPTIVE_BOUND_CANDIDATES")) {
        config.adaptive_candidate_limit = std::max<size_t>(1, static_cast<size_t>(std::stoull(e)));
    }
    config.adaptive_trace_bounds = env_flag_enabled("PYMATCHING_ADAPTIVE_TRACE_BOUNDS");
    if (const char* e = std::getenv("PYMATCHING_ADAPTIVE_BOUND_STRATEGY")) {
        std::string strategy(e);
        if (strategy == "breakpoint" || strategy == "breakpoints") {
            config.adaptive_breakpoint_strategy = true;
        } else if (strategy == "trial" || strategy == "grid") {
            config.adaptive_breakpoint_strategy = false;
        } else {
            throw std::invalid_argument("unknown PYMATCHING_ADAPTIVE_BOUND_STRATEGY: " + strategy);
        }
    }
    if (env_flag_enabled("PYMATCHING_ADAPTIVE_BREAKPOINT_BOUNDS")) {
        config.adaptive_breakpoint_strategy = true;
    }
    config.adaptive_enforce_monotone_bounds = env_flag_enabled("PYMATCHING_ADAPTIVE_ENFORCE_MONOTONE_BOUNDS");
    // Correctness-oriented adaptive defaults: keep the phi certificate and
    // the legacy inter-level gap unless explicitly disabled.  The positive
    // ENFORCE_* variables are still accepted for old scripts, but are now
    // redundant with the default.
    config.adaptive_enforce_legacy_gap_bound =
        !env_flag_enabled("PYMATCHING_ADAPTIVE_DISABLE_LEGACY_GAP_FILTER") ||
        env_flag_enabled("PYMATCHING_ADAPTIVE_ENFORCE_LEGACY_GAP_FILTER");
    config.adaptive_enforce_phi_buffer_certificate =
        !env_flag_enabled("PYMATCHING_ADAPTIVE_DISABLE_PHI_CERTIFICATE") ||
        env_flag_enabled("PYMATCHING_ADAPTIVE_ENFORCE_PHI_CERTIFICATE");
    config.adaptive_min_buffer_ratio_first_level_only =
        env_flag_enabled("PYMATCHING_ADAPTIVE_MIN_BUFFER_RATIO_FIRST_LEVEL_ONLY");
    if (const char* e = std::getenv("PYMATCHING_ADAPTIVE_MIN_BUFFER_RATIO")) {
        config.adaptive_min_buffer_to_diameter_ratio = std::stod(e);
    } else if (!env_flag_enabled("PYMATCHING_ADAPTIVE_DISABLE_MIN_BUFFER_RATIO")) {
        if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_MIN_BUFFER_RATIO")) {
            config.adaptive_min_buffer_to_diameter_ratio = std::stod(e);
        } else {
            config.adaptive_min_buffer_to_diameter_ratio = 2.2;
        }
    }
    if (const char* e = std::getenv("PYMATCHING_ADAPTIVE_PHI_FLOOR")) {
        config.adaptive_phi_floor = std::stod(e);
    } else if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE_FLOOR")) {
        config.adaptive_phi_floor = std::stod(e);
    }
    if (const char* e = std::getenv("PYMATCHING_ADAPTIVE_OBJECTIVE")) {
        std::string objective(e);
        if (objective == "ideal_events" || objective == "level_critical_events" ||
            objective == "parallel_ideal_events" || objective == "levelmax_sum_events" ||
            objective == "critical_events") {
            config.adaptive_objective_minimize_ideal_events = true;
            config.adaptive_objective_minimize_max_events = false;
            config.adaptive_objective_minimize_max_active = false;
        } else if (objective == "max_events" || objective == "min_max_events" ||
            objective == "event_count" || objective == "events") {
            config.adaptive_objective_minimize_ideal_events = false;
            config.adaptive_objective_minimize_max_events = true;
            config.adaptive_objective_minimize_max_active = false;
        } else if (objective == "max_active" || objective == "min_max_active" ||
                   objective == "max_cluster_active") {
            config.adaptive_objective_minimize_ideal_events = false;
            config.adaptive_objective_minimize_max_active = true;
            config.adaptive_objective_minimize_max_events = false;
        } else if (objective == "cluster_count" || objective == "clusters" ||
                   objective == "legacy") {
            config.adaptive_objective_minimize_ideal_events = false;
            config.adaptive_objective_minimize_max_active = false;
            config.adaptive_objective_minimize_max_events = false;
        } else {
            throw std::invalid_argument("unknown PYMATCHING_ADAPTIVE_OBJECTIVE: " + objective);
        }
    }
    if (const char* e = std::getenv("PYMATCHING_ADAPTIVE_LOOKAHEAD_LEVELS")) {
        config.adaptive_lookahead_levels = static_cast<size_t>(std::stoull(e));
    }
    if (const char* e = std::getenv("PYMATCHING_ADAPTIVE_IDEAL_EVENTS_WEIGHT_POWER")) {
        config.adaptive_ideal_events_weight_power = std::stod(e);
    } else if (const char* e = std::getenv("PYMATCHING_ADAPTIVE_LEVEL_EVENT_WEIGHT_POWER")) {
        config.adaptive_ideal_events_weight_power = std::stod(e);
    }
    if (env_flag_enabled("PYMATCHING_ADAPTIVE_LEVEL1_REOPTIMIZE") ||
        env_flag_enabled("PYMATCHING_LEVEL1_REOPTIMIZE_BOUNDS")) {
        config.adaptive_level1_reoptimize_bounds = true;
    }
    if (const char* e = std::getenv("PYMATCHING_ADAPTIVE_LEVEL1_MAX_DIAMETER_RATIO")) {
        config.adaptive_level1_reoptimize_bounds = true;
        config.adaptive_level1_max_diameter_ratio = std::stod(e);
    }
    if (const char* e = std::getenv("PYMATCHING_ADAPTIVE_LEVEL1_MAX_DIAMETER")) {
        config.adaptive_level1_reoptimize_bounds = true;
        config.adaptive_level1_max_diameter = static_cast<pm::cumulative_time_int>(std::stoull(e));
    }
    if (const char* e = std::getenv("PYMATCHING_ADAPTIVE_LEVEL1_CANDIDATES")) {
        config.adaptive_level1_reoptimize_bounds = true;
        config.adaptive_level1_candidate_limit = std::max<size_t>(1, static_cast<size_t>(std::stoull(e)));
    }
    if (env_flag_enabled("PYMATCHING_ADAPTIVE_DEFER_HEAVY_CLUSTERS")) {
        config.adaptive_defer_heavy_clusters = true;
    }
    if (const char* e = std::getenv("PYMATCHING_ADAPTIVE_DEFER_MAX_ACTIVE")) {
        config.adaptive_defer_heavy_clusters = true;
        config.adaptive_defer_max_active = static_cast<size_t>(std::stoull(e));
    }
    if (const char* e = std::getenv("PYMATCHING_ADAPTIVE_DEFER_MAX_EVENTS")) {
        config.adaptive_defer_heavy_clusters = true;
        config.adaptive_defer_max_events = static_cast<uint64_t>(std::stoull(e));
    }
    if (const char* e = std::getenv("PYMATCHING_LEVEL1_SHORT_PAIR_MAX_DISTANCE")) {
        config.level1_short_pair_acceptance = true;
        config.level1_short_pair_max_distance = static_cast<pm::cumulative_time_int>(std::stoull(e));
    }
    if (const char* e = std::getenv("PYMATCHING_LEVEL1_SHORT_PAIR_DISTANCE_RATIO")) {
        config.level1_short_pair_acceptance = true;
        config.level1_short_pair_distance_ratio = std::stod(e);
    }
    if (const char* e = std::getenv("PYMATCHING_LEVEL1_SHORT_PAIR_EXTERNAL_GUARD")) {
        config.level1_short_pair_external_guard = static_cast<pm::cumulative_time_int>(std::stoull(e));
        config.level1_short_pair_external_guard_explicit = true;
    }
    if (env_flag_enabled("PYMATCHING_LEVEL1_SHORT_PAIR_DISABLE_MUTUAL_NEAREST")) {
        config.level1_short_pair_require_mutual_nearest = false;
    }

    // Convenience switch: reuse the existing path-bound schedule options as
    // adaptive d_k,b_k constraints.  This keeps the old knobs available without
    // forcing them on for pure cluster-count maximization experiments.
    if (env_flag_enabled("PYMATCHING_ADAPTIVE_USE_CLUSTER_PATH_CONSTRAINTS")) {
        config.adaptive_enforce_legacy_gap_bound |=
            env_flag_enabled("PYMATCHING_CLUSTER_PATH_LEGACY_GAP_FILTER");
        config.adaptive_min_buffer_ratio_first_level_only |=
            env_flag_enabled("PYMATCHING_CLUSTER_PATH_MIN_BUFFER_RATIO_FIRST_LEVEL_ONLY");
        if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_MIN_BUFFER_RATIO")) {
            config.adaptive_min_buffer_to_diameter_ratio = std::stod(e);
        } else if (!(config.adaptive_min_buffer_to_diameter_ratio > 0.0)) {
            config.adaptive_min_buffer_to_diameter_ratio = 2.2;
        }
        if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_PHI_FLOOR")) {
            config.adaptive_phi_floor = std::stod(e);
        }
    }

    std::cerr << "ADAPTIVE_CLUSTER_BOUNDS_CONFIG"
              << " enabled=1"
              << " candidates=" << config.adaptive_candidate_limit
              << " strategy=" << (config.adaptive_breakpoint_strategy ? "breakpoint" : "trial")
              << " trace=" << (config.adaptive_trace_bounds ? 1 : 0)
              << " enforce_monotone=" << (config.adaptive_enforce_monotone_bounds ? 1 : 0)
              << " enforce_legacy_gap=" << (config.adaptive_enforce_legacy_gap_bound ? 1 : 0)
              << " enforce_phi=" << (config.adaptive_enforce_phi_buffer_certificate ? 1 : 0)
              << " min_buffer_ratio=" << config.adaptive_min_buffer_to_diameter_ratio
              << " min_buffer_ratio_first_level_only=" << (config.adaptive_min_buffer_ratio_first_level_only ? 1 : 0)
              << " phi_floor=" << config.adaptive_phi_floor
              << " objective=" << (config.adaptive_objective_minimize_ideal_events ? "ideal_events" : (config.adaptive_objective_minimize_max_events ? "max_events" : (config.adaptive_objective_minimize_max_active ? "max_active" : "cluster_count")))
              << " lookahead=" << config.adaptive_lookahead_levels
              << " ideal_events_weight_power=" << config.adaptive_ideal_events_weight_power
              << " level1_reopt=" << (config.adaptive_level1_reoptimize_bounds ? 1 : 0)
              << " level1_reopt_max_ratio=" << config.adaptive_level1_max_diameter_ratio
              << " level1_reopt_max_d=" << config.adaptive_level1_max_diameter
              << " level1_reopt_candidates=" << config.adaptive_level1_candidate_limit
              << " defer_heavy=" << (config.adaptive_defer_heavy_clusters ? 1 : 0)
              << " defer_max_active=" << config.adaptive_defer_max_active
              << " defer_max_events=" << config.adaptive_defer_max_events
              << " level1_short_pair=" << (config.level1_short_pair_acceptance ? 1 : 0)
              << " level1_short_pair_max_distance=" << config.level1_short_pair_max_distance
              << " level1_short_pair_distance_ratio=" << config.level1_short_pair_distance_ratio
              << " level1_short_pair_external_guard=" << config.level1_short_pair_external_guard
              << " level1_short_pair_external_guard_explicit=" << (config.level1_short_pair_external_guard_explicit ? 1 : 0)
              << " level1_short_pair_mutual_nearest=" << (config.level1_short_pair_require_mutual_nearest ? 1 : 0)
              << "\n";
}

void run_and_print_growing_only_benchmark(
    const std::string& case_name,
    pm::UserGraph& graph,
    const std::vector<Shot>& shots,
    int reps,
    const pm::ProcessingClusterConfig& cluster_config,
    size_t max_levels_to_run) {
    const auto& matching_graph = graph.get_matching_graph_for_parallel_clustering();
    auto* cache = graph.get_processing_cluster_graph_cache();
    pm::ensure_processing_cluster_radius_neighbors_precomputed(matching_graph, *cache, cluster_config);
    size_t effective_max_level = pm::processing_cluster_effective_max_level(cluster_config);
    if (max_levels_to_run == 0 || max_levels_to_run > effective_max_level) {
        max_levels_to_run = effective_max_level;
    }
    std::vector<GrowingOnlyTotals> totals(max_levels_to_run + 1);
    std::vector<pm::ProcessingClusterBounds> bounds(max_levels_to_run + 1);
    for (size_t level = 1; level <= max_levels_to_run; level++) {
        bounds[level] = pm::processing_cluster_bounds_for_level(level, cluster_config);
    }

    const bool dump_examples = env_flag_enabled("PYMATCHING_GROWING_ONLY_DUMP_CLUSTERS");
    size_t dump_shot = 0;
    if (const char* e = std::getenv("PYMATCHING_GROWING_ONLY_DUMP_SHOT")) {
        dump_shot = static_cast<size_t>(std::stoull(e));
    }
    size_t max_dump_clusters = 8;
    if (const char* e = std::getenv("PYMATCHING_GROWING_ONLY_DUMP_MAX_CLUSTERS")) {
        max_dump_clusters = static_cast<size_t>(std::stoull(e));
    }

    for (int r = 0; r < reps; r++) {
        for (const auto& shot : shots) {
            size_t shot_index = static_cast<size_t>(&shot - shots.data());
            std::vector<uint64_t> residual = shot.hits;
            pm::ProcessingClusterConfig shot_cluster_config = cluster_config;
            for (size_t level = 1; level <= max_levels_to_run; level++) {
                auto& t = totals[level];
                t.input_active += residual.size();
                if (residual.empty()) {
                    continue;
                }
                auto start = eventcount_optional_now();
                auto result = pm::run_growing_only_clustering_with_adaptive_bounds(
                    matching_graph, std::move(residual), level, shot_cluster_config, cache);
                bounds[level] = result.bounds;
                uint64_t wall = ns_since(start);
                t.wall_ns += wall;
                t.event_generation_ns += result.stats.event_generation_wall_ns;
                t.event_processing_ns += result.stats.event_processing_wall_ns;
                t.generated_events += result.stats.generated_pair_collision_events;
                t.processed_events += result.stats.processed_pair_collision_events;
                t.internal_collisions += result.stats.internal_collisions_recorded;
                t.external_collisions += result.stats.external_collisions_recorded;
                t.union_count += result.stats.union_count;
                t.accepted_clusters += result.accepted_clusters.size();
                t.residual_active += result.residual_active_detector_count();
                t.removed_or_hidden_clusters += result.stats.removed_or_hidden_clusters;
                t.parity_rejected_clusters += result.stats.parity_rejected_clusters;
                t.diameter_rejected_clusters += result.stats.diameter_rejected_clusters;
                t.max_active_rejected_clusters += result.stats.max_active_rejected_clusters;
                t.forced_by_max_level_clusters += result.stats.forced_by_max_level_clusters;
                t.sparse_frontier_events += result.stats.sparse_frontier_events;
                t.processed_frontier_events += result.stats.processed_frontier_events;
                t.generated_collision_events += result.stats.generated_collision_events;
                t.processed_collision_events += result.stats.processed_collision_events;
                t.exact_internal_pair_checks += result.stats.exact_internal_pair_checks;
                t.exact_external_pair_checks += result.stats.exact_external_pair_checks;
                t.exact_diameter_checks += result.stats.exact_diameter_checks;
                t.sparse_frontier_used += result.stats.used_sparse_frontier ? 1 : 0;
                t.component_lookup_used += result.stats.used_component_lookup ? 1 : 0;
                t.active_pair_table_used += result.stats.used_active_pair_table ? 1 : 0;
                for (const auto& c : result.accepted_clusters) {
                    t.accepted_active += c.active_detectors.size();
                    t.accepted_cluster_size_sum += c.active_detectors.size();
                    t.accepted_cluster_max_size = std::max<uint64_t>(t.accepted_cluster_max_size, c.active_detectors.size());
                    t.max_diameter = std::max(t.max_diameter, c.diameter);
                    t.max_buffer = std::max(t.max_buffer, c.buffer);
                    t.max_ready_time = std::max(t.max_ready_time, c.ready_time);
                    if (c.first_external_collision_time < pm::GROWING_ONLY_INF_DISTANCE) {
                        t.max_first_external_collision_time = std::max(
                            t.max_first_external_collision_time, c.first_external_collision_time);
                    }
                }
                if (dump_examples && shot_index == dump_shot && t.examples_printed < max_dump_clusters) {
                    for (const auto& c : result.accepted_clusters) {
                        if (t.examples_printed >= max_dump_clusters) {
                            break;
                        }
                        std::cerr << "GROWING_ONLY_CLUSTER_EXAMPLE case=" << case_name
                                  << " shot=" << shot_index
                                  << " level=" << level
                                  << " d=" << result.bounds.diameter_bound
                                  << " b=" << result.bounds.buffer_bound
                                  << " size=" << c.active_detectors.size()
                                  << " diameter=" << c.diameter
                                  << " buffer=" << c.buffer
                                  << " ready_time=" << c.ready_time
                                  << " max_internal_pair_collision_time=" << c.max_internal_pair_collision_time
                                  << " max_internal_pair_distance=" << c.max_internal_pair_distance
                                  << " first_external_collision_time=" << c.first_external_collision_time
                                  << " nearest_boundary_match_distance=" << c.nearest_boundary_match_distance
                                  << " parity_can_stop=" << (c.parity_can_stop_locally ? 1 : 0)
                                  << " forced_by_max_level=" << (c.forced_by_max_level ? 1 : 0)
                                  << " internal_collisions=" << c.internal_collisions_recorded
                                  << " external_collisions=" << c.external_collisions_recorded
                                  << " detectors=";
                        for (size_t k = 0; k < c.active_detectors.size(); k++) {
                            if (k) std::cerr << ":";
                            std::cerr << c.active_detectors[k];
                        }
                        std::cerr << "\n";
                        t.examples_printed++;
                    }
                }
                residual = result.residual_active_detectors;
            }
        }
    }

    double shot_reps = static_cast<double>(shots.size()) * static_cast<double>(reps);
    for (size_t level = 1; level <= max_levels_to_run; level++) {
        const auto& t = totals[level];
        const auto& b = bounds[level];
        const double accepted_clusters = shot_reps ? static_cast<double>(t.accepted_clusters) / shot_reps : 0.0;
        const double avg_cluster_size = t.accepted_clusters ? static_cast<double>(t.accepted_cluster_size_sum) / static_cast<double>(t.accepted_clusters) : 0.0;
        const double removed_fraction = shots.empty() ? 0.0 : static_cast<double>(t.accepted_active) / static_cast<double>(std::max<uint64_t>(1, t.accepted_active + t.residual_active));
        std::cerr << std::fixed << std::setprecision(6)
                  << "GROWING_ONLY_SUMMARY case=" << case_name
                  << " level=" << level
                  << " d=" << b.diameter_bound
                  << " b=" << b.buffer_bound
                  << " shots=" << static_cast<uint64_t>(shot_reps)
                  << " wall_ms_per_shot=" << static_cast<double>(t.wall_ns) / shot_reps / 1e6
                  << " event_generation_ms_per_shot=" << static_cast<double>(t.event_generation_ns) / shot_reps / 1e6
                  << " event_processing_ms_per_shot=" << static_cast<double>(t.event_processing_ns) / shot_reps / 1e6
                  << " accepted_clusters_per_shot=" << accepted_clusters
                  << " input_active_per_shot=" << static_cast<double>(t.input_active) / shot_reps
                  << " residual_active_per_shot=" << static_cast<double>(t.residual_active) / shot_reps
                  << " accepted_active_per_shot=" << static_cast<double>(t.accepted_active) / shot_reps
                  << " removed_active_fraction=" << removed_fraction
                  << " avg_cluster_size=" << avg_cluster_size
                  << " max_cluster_size=" << t.accepted_cluster_max_size
                  << " max_diameter=" << t.max_diameter
                  << " max_buffer=" << t.max_buffer
                  << " max_ready_time=" << t.max_ready_time
                  << " max_first_external_collision_time=" << t.max_first_external_collision_time
                  << " internal_collisions_per_shot=" << static_cast<double>(t.internal_collisions) / shot_reps
                  << " external_collisions_per_shot=" << static_cast<double>(t.external_collisions) / shot_reps
                  << " generated_events_per_shot=" << static_cast<double>(t.generated_events) / shot_reps
                  << " processed_events_per_shot=" << static_cast<double>(t.processed_events) / shot_reps
                  << " unions_per_shot=" << static_cast<double>(t.union_count) / shot_reps
                  << " removed_or_hidden_clusters_per_shot=" << static_cast<double>(t.removed_or_hidden_clusters) / shot_reps
                  << " parity_rejected_clusters_per_shot=" << static_cast<double>(t.parity_rejected_clusters) / shot_reps
                  << " diameter_rejected_clusters_per_shot=" << static_cast<double>(t.diameter_rejected_clusters) / shot_reps
                  << " max_active_rejected_clusters_per_shot=" << static_cast<double>(t.max_active_rejected_clusters) / shot_reps
                  << " forced_by_max_level_clusters_per_shot=" << static_cast<double>(t.forced_by_max_level_clusters) / shot_reps
                  << " sparse_frontier_used_per_shot=" << static_cast<double>(t.sparse_frontier_used) / shot_reps
                  << " component_lookup_used_per_shot=" << static_cast<double>(t.component_lookup_used) / shot_reps
                  << " active_pair_table_used_per_shot=" << static_cast<double>(t.active_pair_table_used) / shot_reps
                  << " sparse_frontier_events_per_shot=" << static_cast<double>(t.sparse_frontier_events) / shot_reps
                  << " processed_frontier_events_per_shot=" << static_cast<double>(t.processed_frontier_events) / shot_reps
                  << " generated_collision_events_per_shot=" << static_cast<double>(t.generated_collision_events) / shot_reps
                  << " processed_collision_events_per_shot=" << static_cast<double>(t.processed_collision_events) / shot_reps
                  << " exact_internal_pair_checks_per_shot=" << static_cast<double>(t.exact_internal_pair_checks) / shot_reps
                  << " exact_external_pair_checks_per_shot=" << static_cast<double>(t.exact_external_pair_checks) / shot_reps
                  << " exact_diameter_checks_per_shot=" << static_cast<double>(t.exact_diameter_checks) / shot_reps
                  << "\n";
    }
}

Totals run_parallel(
    pm::UserGraph& graph,
    const std::vector<Shot>& shots,
    int reps,
    size_t workers,
    size_t ideal_workers,
    double beta,
    double gamma,
    double lambda,
    bool enable_checkpoints,
    const std::vector<uint64_t>* global_observable_predictions = nullptr,
    const std::vector<pm::ProcessingClusterBounds>* explicit_bounds = nullptr,
    const std::vector<size_t>* max_active_detector_caps = nullptr,
    const std::string* case_name_for_cluster_csv = nullptr,
    const CaseData* case_data_for_cluster_csv = nullptr) {
    (void)graph.get_processing_cluster_graph_cache();
    Totals t;
    pm::LockstepSchedulerConfig config;
    config.cluster_config.beta = beta;
    config.cluster_config.gamma = gamma;
    config.cluster_config.lambda = lambda;
    config.cluster_config.max_level = explicit_bounds == nullptr ? 64 : explicit_bounds->size();
    if (explicit_bounds != nullptr) {
        config.cluster_config.explicit_bounds_by_level = *explicit_bounds;
        config.cluster_config.enforce_stopping_lemma_bounds =
            env_flag_enabled("PYMATCHING_ENFORCE_STOPPING_LEMMA_BOUNDS");
    }
    if (max_active_detector_caps != nullptr) {
        config.cluster_config.max_active_detectors_by_level = *max_active_detector_caps;
    }
    apply_adaptive_cluster_bounds_env(config.cluster_config);
    std::map<std::vector<uint64_t>, uint64_t> adaptive_event_count_cache;
    std::map<std::pair<std::vector<uint64_t>, pm::cumulative_time_int>, pm::AdaptiveEventTiming> adaptive_event_timing_cache;
    if (config.cluster_config.adaptive_objective_minimize_max_events ||
        config.cluster_config.adaptive_objective_minimize_ideal_events) {
        auto& adaptive_mwpm = graph.get_mwpm();
        config.cluster_config.adaptive_event_count_callback = [&](const std::vector<uint64_t>& detectors) -> uint64_t {
            std::vector<uint64_t> key = detectors;
            std::sort(key.begin(), key.end());
            key.erase(std::unique(key.begin(), key.end()), key.end());
            auto it = adaptive_event_count_cache.find(key);
            if (it != adaptive_event_count_cache.end()) {
                return it->second;
            }
            uint64_t events = count_mwpm_events_only(adaptive_mwpm, key);
            adaptive_event_count_cache.emplace(std::move(key), events);
            return events;
        };
        config.cluster_config.adaptive_event_timing_callback = [&](
            const std::vector<uint64_t>& detectors,
            pm::cumulative_time_int cutoff_algorithmic_time) -> pm::AdaptiveEventTiming {
            std::vector<uint64_t> key = detectors;
            std::sort(key.begin(), key.end());
            key.erase(std::unique(key.begin(), key.end()), key.end());
            auto cache_key = std::make_pair(key, cutoff_algorithmic_time);
            auto it = adaptive_event_timing_cache.find(cache_key);
            if (it != adaptive_event_timing_cache.end()) {
                return it->second;
            }
            pm::AdaptiveEventTiming timing = count_mwpm_event_timing(adaptive_mwpm, key, cutoff_algorithmic_time);
            adaptive_event_timing_cache.emplace(std::move(cache_key), timing);
            adaptive_event_count_cache.emplace(key, timing.raw_event_count);
            return timing;
        };
    }
    config.allow_global_fallback = false;
    config.num_workers = workers;
    config.ideal_cluster_worker_count = ideal_workers;
    (void)enable_checkpoints;
    config.enable_single_root_direct_decode_fast_path = false;
    config.enable_root_direct_decode_fast_path = false;
    config.enable_root_global_decode_fast_path = false;

    // Current policy path only:
    //   full-graph workers
    //   + growing-only processing clusters
    //   + parentless broadcast import
    //   + level-by-level pipeline
    //
    // Historical cut-out subgraph, root/global bypass, per-child immediate import,
    // and overlap-checkpoint experiment paths were intentionally removed from the
    // benchmark driver so that one source path corresponds to the reported policy.
    config.enable_full_graph_worker_subgraphs = true;
    config.exclude_full_graph_worker_distribution_from_runtime = true;
    config.enable_level_batched_scheduler = true;
    config.enable_parentless_broadcast_import = true;
    config.enable_growing_only_processing_clusters = true;
    config.enable_parentless_level_pipeline = true;
    config.cluster_config.root_direct_min_level = 1;
    config.cluster_config.root_direct_skip_diameter_check = false;
    config.enable_overlap_checkpoints = false;

    auto warmup_start = eventcount_optional_now();

    // One-time graph/parameter preprocessing.  Radius-neighbor lists are built
    // for levels 1..inferred_root_level before the per-shot timer starts.
    auto* precompute_cache = graph.get_processing_cluster_graph_cache();
    const auto& precompute_matching_graph = graph.get_matching_graph_for_parallel_clustering();
    pm::ensure_processing_cluster_radius_neighbors_precomputed(
        precompute_matching_graph, *precompute_cache, config.cluster_config);
    if (workers > 0) {
        pm::prewarm_persistent_full_graph_worker_states(graph, workers);
    }
    if (env_flag_enabled("PYMATCHING_PARENTLESS_ACTUAL_PARALLEL_LEVEL") && workers > 1) {
        pm::prewarm_parentless_level_worker_pool(workers);
    }

    int warmup_reps = 0;
    if (const char* e = std::getenv("PYMATCHING_WARMUP_REPS")) {
        warmup_reps = std::max(0, std::stoi(e));
    }
    size_t required_ready_states = workers;
    if (const char* e = std::getenv("PYMATCHING_REQUIRED_READY_STATES")) {
        required_ready_states = static_cast<size_t>(std::stoull(e));
    }
    const bool measure_bulk_execution_wall = env_flag_enabled("PYMATCHING_MEASURE_EXECUTION_WALL");

    for (int wr = 0; wr < warmup_reps; wr++) {
        for (const auto& warmup_shot : shots) {
            try {
                (void)pm::lockstep_hierarchical_decode(graph, warmup_shot.hits, config);
            } catch (const std::exception& ex) {
                std::cerr << "WARMUP_EXCEPTION what=" << ex.what() << "\n";
                throw;
            }
        }
    }
    if (workers > 0) {
        pm::wait_for_persistent_full_graph_worker_states_ready(graph, required_ready_states);
    }
    t.warmup_wall_ns = ns_since(warmup_start);
    std::cerr << "WARMUP_DONE reps=" << warmup_reps
              << " required_ready_states=" << required_ready_states
              << " wall_ms=" << (static_cast<double>(t.warmup_wall_ns) / 1e6) << "\n";

    auto execution_start = eventcount_optional_now();

    for (int r = 0; r < reps; r++) {
        for (const auto& shot : shots) {
            size_t shot_index = static_cast<size_t>(&shot - shots.data());
            try {
                bool dumped_hierarchy = false;
                bool dump_hierarchy_all = false;
                if (const char* dump_all_env = std::getenv("PYMATCHING_DUMP_HIERARCHY_ALL_SHOTS")) {
                    dump_hierarchy_all = std::string(dump_all_env) == "1";
                }
                bool dump_hierarchy_this_shot = dump_hierarchy_all;
                if (const char* dump_env = std::getenv("PYMATCHING_DUMP_HIERARCHY_SHOT")) {
                    size_t dump_idx = static_cast<size_t>(std::stoull(dump_env));
                    size_t shot_index_for_dump = static_cast<size_t>(&shot - shots.data());
                    dump_hierarchy_this_shot = dump_hierarchy_this_shot || shot_index_for_dump == dump_idx;
                }
                if (dump_hierarchy_this_shot) {
                    size_t shot_index_for_dump = static_cast<size_t>(&shot - shots.data());
                    pm::ProcessingClusterProfilingStats dump_stats;
                    auto* cache = graph.get_processing_cluster_graph_cache();
                    const auto& matching_graph = graph.get_matching_graph_for_parallel_clustering();
                    pm::ensure_processing_cluster_radius_neighbors_precomputed(matching_graph, *cache, config.cluster_config);
                    auto h = pm::build_processing_cluster_hierarchy(
                        matching_graph, shot.hits, config.cluster_config, cache, &dump_stats);
                    std::cerr << "HIERARCHY shot=" << shot_index_for_dump << " clusters=" << h.clusters.size() << "\n";
                    for (size_t level = 1; level < h.cluster_ids_by_level.size(); level++) {
                        if (h.cluster_ids_by_level[level].empty()) {
                            continue;
                        }
                        size_t active_total = 0;
                        size_t max_active = 0;
                        pm::cumulative_time_int max_diameter_seen = 0;
                        pm::cumulative_time_int level_diameter_bound = 0;
                        pm::cumulative_time_int level_buffer_bound = 0;
                        for (auto cluster_id : h.cluster_ids_by_level[level]) {
                            const auto& c = h.clusters[cluster_id];
                            level_diameter_bound = std::max(level_diameter_bound, c.diameter_bound);
                            level_buffer_bound = std::max(level_buffer_bound, c.buffer_bound);
                            active_total += c.active_detectors.size();
                            max_active = std::max(max_active, c.active_detectors.size());
                            max_diameter_seen = std::max(max_diameter_seen, c.diameter);
                        }
                        std::cerr << "HIERARCHY_LEVEL shot=" << shot_index_for_dump
                                  << " level=" << level
                                  << " d=" << level_diameter_bound
                                  << " b=" << level_buffer_bound
                                  << " clusters=" << h.cluster_ids_by_level[level].size()
                                  << " active_total=" << active_total
                                  << " max_active=" << max_active
                                  << " max_diameter=" << max_diameter_seen << "\n";
                    }
                    bool dump_clusters = !dump_hierarchy_all;
                    if (const char* dump_cl_env = std::getenv("PYMATCHING_DUMP_HIERARCHY_CLUSTERS")) {
                        dump_clusters = std::string(dump_cl_env) == "1";
                    }
                    if (dump_clusters) {
                        for (const auto& c : h.clusters) {
                            bool has1192 = std::find(c.influence_vertices.begin(), c.influence_vertices.end(), (size_t)1192) != c.influence_vertices.end();
                            bool active1192 = std::find(c.active_detectors.begin(), c.active_detectors.end(), (uint64_t)1192) != c.active_detectors.end();
                            std::cerr << "CL shot=" << shot_index_for_dump << " id=" << c.id << " level=" << c.level << " parent=" << c.parent_id
                                      << " active=" << c.active_detectors.size()
                                      << " influence=" << c.influence_vertices.size()
                                      << " boundary=" << c.boundary_endpoint_vertices.size()
                                      << " has1192=" << has1192 << " active1192=" << active1192
                                      << " full=" << c.influence_is_full_graph << "\n";
                        }
                    }
                    dumped_hierarchy = true;
                }
                if (dumped_hierarchy) {
                    if (const char* only_dump_env = std::getenv("PYMATCHING_ONLY_DUMP_HIERARCHY")) {
                        if (std::string(only_dump_env) == "1") {
                            continue;
                        }
                    }
                }
                const bool debug_trace_shots = env_flag_enabled("PYMATCHING_DEBUG_TRACE_SHOTS");
                if (debug_trace_shots) {
                    std::cerr << "PYMATCHING DEBUG SHOT START rep=" << r
                              << " shot=" << shot_index
                              << " active=" << shot.hits.size() << std::endl;
                }
                auto start = eventcount_optional_now();
                auto result = pm::lockstep_hierarchical_decode(graph, shot.hits, config);
                auto shot_decode_ns = ns_since(start);
                if (debug_trace_shots) {
                    std::cerr << "PYMATCHING DEBUG SHOT DONE rep=" << r
                              << " shot=" << shot_index
                              << " active=" << shot.hits.size()
                              << " wall_ms=" << (static_cast<double>(shot_decode_ns) / 1e6)
                              << " clusters=" << result.hierarchy.clusters.size()
                              << " max_level=" << result.profiling_stats.max_level << std::endl;
                }
                if (case_name_for_cluster_csv != nullptr && case_data_for_cluster_csv != nullptr) {
                    if (const char* cluster_csv_env = std::getenv("PYMATCHING_EVENTCOUNT_CLUSTER_CSV")) {
                        if (*cluster_csv_env != '\0') {
                            append_cluster_event_csv_rows(
                                cluster_csv_env, *case_name_for_cluster_csv, *case_data_for_cluster_csv, r, shot_index, result);
                        }
                    }
                    append_fixed_cluster_detail_csv_rows_if_requested(
                        *case_name_for_cluster_csv, *case_data_for_cluster_csv, r, shot_index, result);
                }
                if (!measure_bulk_execution_wall) {
                    t.wall_ns += shot_decode_ns;
                }
                auto& s = result.profiling_stats;
                t.clustering_ns += s.clustering_wall_ns;
                t.component_ns += s.clustering_component_construction_wall_ns;
                t.diameter_ns += s.clustering_diameter_check_wall_ns;
                t.precomputable_distance_lookup_ns += s.clustering_precomputable_distance_lookup_wall_ns;
                t.influence_ns += s.clustering_influence_region_wall_ns;
                t.parent_ns += s.clustering_parent_assignment_wall_ns;
                t.used_active_pair_table += s.growing_only_used_active_pair_table;
                t.used_sparse_frontier += s.growing_only_used_sparse_frontier;
                t.used_component_lookup += s.growing_only_used_component_lookup;
                t.forced_max_level_shortcut_count += s.growing_only_forced_max_level_shortcut_count;
                t.processed_frontier_events += s.growing_only_processed_frontier_events;
                t.generated_collision_events += s.growing_only_generated_collision_events;
                t.processed_collision_events += s.growing_only_processed_collision_events;
                t.exact_internal_pair_checks += s.growing_only_exact_internal_pair_checks;
                t.exact_external_pair_checks += s.growing_only_exact_external_pair_checks;
                t.exact_diameter_checks += s.growing_only_exact_diameter_checks;
                t.max_component_active_size = std::max<uint64_t>(t.max_component_active_size, s.growing_only_max_component_active_size);
                t.max_residual_active_size = std::max<uint64_t>(t.max_residual_active_size, s.growing_only_max_residual_active_size);
                t.subgraph_build_ns += s.subgraph_build_wall_ns;
                t.worker_graph_distribution_excluded_ns += s.worker_graph_distribution_wall_ns;
                t.initial_setup_ns += s.initial_mwpm_setup_wall_ns;
                t.scheduler_ns += s.scheduler_wall_ns;
                t.child_import_ns += s.child_import_wall_ns;
                t.checkpoint_capture_ns += s.overlap_checkpoint_capture_wall_ns;
                t.checkpoint_restore_ns += s.overlap_checkpoint_restore_wall_ns;
                t.cluster_step_ns += s.cluster_step_wall_ns;
                t.max_level_cluster_ns += max_level_cluster_work_ns(s);
                t.max_level_lifecycle_ns += max_level_cluster_lifecycle_ns(s);
                t.parallel_policy_critical_path_lifecycle_ns += parallel_policy_critical_path_lifecycle_ns(s);
                t.ideal_parallel_time_ns += ideal_parallel_time_ns(s);
                t.ideal_runtime_preprocessing_plus_parallel_time_ns +=
                    ideal_runtime_preprocessing_plus_parallel_time_ns(s);
                t.runtime_preprocessing_plus_parallel_sparse_blossom_ns += runtime_preprocessing_plus_parallel_sparse_blossom_ns(s);
                t.runtime_preprocessing_excluding_precompute_ns += runtime_preprocessing_excluding_precompute_ns(s);
                t.runtime_parallel_sparse_blossom_ns += parallel_policy_critical_path_lifecycle_ns(s);
                t.max_level_lifecycle_initial_setup_ns += s.max_level_cluster_initial_setup_wall_ns;
                t.max_level_lifecycle_child_import_ns += s.max_level_cluster_child_import_wall_ns;
                t.max_level_lifecycle_checkpoint_capture_ns += s.max_level_cluster_checkpoint_capture_wall_ns;
                t.max_level_lifecycle_checkpoint_restore_ns += s.max_level_cluster_checkpoint_restore_wall_ns;
                t.max_level_lifecycle_mark_stop_ns += s.max_level_cluster_mark_stop_wall_ns;
                t.max_level_lifecycle_cluster_id += s.max_level_cluster_id;
                t.overhead_ns += parallel_overhead_ns(s);
                t.ideal_ns += ideal_parallel_ns(s);
                t.overhead_excluding_precompute_ns += overhead_excluding_precompute_ns(s);
                t.ideal_excluding_precompute_ns += ideal_parallel_excluding_precompute_ns(s);
                t.single_core_preprocessing_excluding_precompute_ns += single_core_preprocessing_excluding_precompute_ns(s);
                t.ideal_worker_cluster_makespan_ns += s.ideal_worker_cluster_makespan_wall_ns;
                t.ideal_worker_communication_estimate_ns += s.ideal_worker_communication_estimate_ns;
                t.multi_worker_ideal_excluding_precompute_with_comm_ns +=
                    multi_worker_ideal_excluding_precompute_with_comm_ns(s);
                t.ideal_worker_count_used += s.ideal_cluster_worker_count_used;
                t.direct_worker_init_ns += s.direct_worker_init_wall_ns;
                t.direct_worker_advance_ns += s.direct_worker_advance_wall_ns;
                t.direct_worker_extract_ns += s.direct_worker_extract_wall_ns;
                t.direct_worker_release_ns += s.direct_worker_release_wall_ns;
                t.direct_worker_total_ns += s.direct_worker_total_wall_ns;
                t.max_direct_worker_init_ns += s.max_direct_worker_init_wall_ns;
                t.max_direct_worker_advance_ns += s.max_direct_worker_advance_wall_ns;
                t.max_direct_worker_extract_ns += s.max_direct_worker_extract_wall_ns;
                t.max_direct_worker_release_ns += s.max_direct_worker_release_wall_ns;
                t.max_direct_worker_total_ns += s.max_direct_worker_total_wall_ns;
                t.direct_worker_active_detectors += s.direct_worker_active_detectors;
                t.max_direct_worker_active_detectors += s.max_direct_worker_active_detectors;
                t.parallel_cluster_events += s.parallel_cluster_events;
                t.parallel_cluster_events_square_sum += static_cast<long double>(s.parallel_cluster_events) * static_cast<long double>(s.parallel_cluster_events);
                t.parallel_level_max_event_sum += s.parallel_level_max_event_sum;
                t.parallel_level_max_event_square_sum += static_cast<long double>(s.parallel_level_max_event_sum) * static_cast<long double>(s.parallel_level_max_event_sum);
                t.parallel_max_cluster_events = std::max<uint64_t>(t.parallel_max_cluster_events, s.parallel_max_cluster_events);
                t.parallel_nonempty_event_levels += s.parallel_nonempty_event_levels;
                const uint64_t parallel_algorithmic_time_sample = static_cast<uint64_t>(
                    std::max<pm::cumulative_time_int>(0, s.parallel_level_critical_path_algorithmic_time));
                t.parallel_critical_path_algorithmic_time += parallel_algorithmic_time_sample;
                t.parallel_critical_path_algorithmic_time_square_sum += static_cast<long double>(parallel_algorithmic_time_sample) * static_cast<long double>(parallel_algorithmic_time_sample);
                t.child_import_count += s.child_import_count;
                t.clusters += s.cluster_count;
                t.max_level += s.max_level;
                t.peak_width += s.max_clusters_in_level;
                t.max_cluster_size += s.max_cluster_active_detectors;
                t.scheduler_batches += s.scheduler_batch_count;
                t.scheduler_cluster_steps += s.scheduler_cluster_step_count;
                uint64_t pred = 0;
                for (size_t k = 0; k < result.root_aggregate_result.obs_crossed.size() && k < 64; k++) {
                    pred ^= (uint64_t)result.root_aggregate_result.obs_crossed[k] << k;
                }
                uint64_t expected_pred_for_log = shot.obs_mask;
                if (global_observable_predictions != nullptr && !global_observable_predictions->empty()) {
                    // Mismatch policy: compare only predicted observable flips.
                    // The full MWPM matching edge set is allowed to differ.
                    expected_pred_for_log = (*global_observable_predictions)[shot_index];
                }
                if (pred != expected_pred_for_log) {
                    t.mistakes++;
                    if (case_name_for_cluster_csv != nullptr && case_data_for_cluster_csv != nullptr) {
                        append_mistake_log_rows_if_requested(
                            *case_name_for_cluster_csv, *case_data_for_cluster_csv, r, shot_index,
                            method_label_from_env("fixed"), pred, expected_pred_for_log, result);
                    }
                }
            } catch (const std::exception& ex) {
                std::cerr << "SHOT_EXCEPTION " << shot_index << "\n";
                std::cerr << "PARALLEL_EXCEPTION shot=" << shot_index
                          << " beta=" << beta << " gamma=" << gamma << " lambda=" << lambda
                          << " hits=" << shot.hits.size()
                          << " what=" << ex.what() << "\n";
                t.exceptions++;
            }
        }
    }
    t.execution_wall_ns = ns_since(execution_start);
    if (measure_bulk_execution_wall) {
        t.wall_ns = t.execution_wall_ns;
    }
    std::cerr << "RUN_DONE bulk_execution_wall=" << (measure_bulk_execution_wall ? 1 : 0)
              << " reps=" << reps
              << " shots=" << shots.size()
              << "\n";
    return t;
}

double standard_error_of_mean(long double sum, long double square_sum, double sample_count) {
    if (sample_count <= 1.0) {
        return 0.0;
    }
    const long double n = static_cast<long double>(sample_count);
    const long double mean = sum / n;
    long double sample_variance = (square_sum - n * mean * mean) / (n - 1.0L);
    if (sample_variance < 0.0L && sample_variance > -1e-9L) {
        sample_variance = 0.0L;
    }
    if (sample_variance <= 0.0L) {
        return 0.0;
    }
    return std::sqrt(static_cast<double>(sample_variance / n));
}

void print_parallel_row(const std::string& case_name, const std::string& mode, const Totals& t, double shots) {
    (void)mode;
    std::cout << std::fixed << std::setprecision(6)
        << case_name << "," << mode
        << "," << (double)t.global_mwpm_events / shots
        << "," << standard_error_of_mean(static_cast<long double>(t.global_mwpm_events), t.global_mwpm_events_square_sum, shots)
        << "," << (double)t.parallel_cluster_events / shots
        << "," << standard_error_of_mean(static_cast<long double>(t.parallel_cluster_events), t.parallel_cluster_events_square_sum, shots)
        << "," << (double)t.parallel_level_max_event_sum / shots
        << "," << standard_error_of_mean(static_cast<long double>(t.parallel_level_max_event_sum), t.parallel_level_max_event_square_sum, shots)
        << "," << t.parallel_max_cluster_events
        << "," << (double)t.parallel_nonempty_event_levels / shots
        << "," << (double)t.parallel_critical_path_algorithmic_time / shots
        << "," << standard_error_of_mean(static_cast<long double>(t.parallel_critical_path_algorithmic_time), t.parallel_critical_path_algorithmic_time_square_sum, shots)
        << "," << (double)t.child_import_count / shots
        << "," << (double)t.clusters / shots
        << "," << (double)t.max_level / shots
        << "," << (double)t.peak_width / shots
        << "," << (double)t.max_cluster_size / shots
        << "," << t.mistakes
        << "," << t.exceptions
        << "\n";
}

void segv_handler(int sig) {
    void* array[80];
    int n = backtrace(array, 80);
    const char msg[] = "\nSIGSEGV backtrace:\n";
    write(2, msg, sizeof(msg)-1);
    backtrace_symbols_fd(array, n, 2);
    _exit(128 + sig);
}

int main(int argc, char** argv) {
    signal(SIGSEGV, segv_handler);
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " reps case_path...\n";
        return 2;
    }
    int reps = std::stoi(argv[1]);
    std::cout << "case,mode,global_mwpm_events_per_shot,global_mwpm_events_std_error_per_shot,parallel_cluster_events_total_per_shot,parallel_cluster_events_total_std_error_per_shot,parallel_level_critical_events_per_shot,parallel_level_critical_events_std_error_per_shot,parallel_critical_events_seen,parallel_nonempty_event_levels_per_shot,parallel_critical_path_algorithmic_time_per_shot,parallel_critical_path_algorithmic_time_std_error_per_shot,child_import_messages_per_shot,avg_clusters,avg_max_level,avg_peak_width,avg_max_cluster_size,mistakes,exceptions\n";
    for (int a = 2; a < argc; a++) {
        std::string path = argv[a];
        auto data = read_case(path);
        // Keep the full case shots available for match-path based hierarchy selection.
        // Low-p smoke tests often set PYMATCHING_BENCH_MAX_SHOTS=1, and the first
        // shot can easily have no detection events.  The benchmark should still be
        // able to infer a stable cluster path unit from the full generated case.
        const auto path_stats_shots = data.shots;
        // Optional debugging aid: run a suffix of the shot list.  This is useful
        // when PYMATCHING_BENCH_MAX_SHOTS=N shows that shot N-1 is pathological;
        // set PYMATCHING_BENCH_SHOT_START=N-1 and PYMATCHING_BENCH_MAX_SHOTS=1
        // to isolate that original shot while still using the full case for
        // match-path based hierarchy parameter selection.
        if (const char* e = std::getenv("PYMATCHING_BENCH_SHOT_START")) {
            size_t shot_start = static_cast<size_t>(std::stoull(e));
            if (shot_start >= data.shots.size()) {
                data.shots.clear();
            } else if (shot_start > 0) {
                data.shots.erase(data.shots.begin(), data.shots.begin() + static_cast<std::ptrdiff_t>(shot_start));
            }
        }
        if (const char* e = std::getenv("PYMATCHING_BENCH_MAX_SHOTS")) {
            size_t max_shots = static_cast<size_t>(std::stoull(e));
            if (max_shots < data.shots.size()) {
                data.shots.resize(max_shots);
            }
        }
        std::string case_name = "r" + std::to_string(data.rounds) + "_d" + std::to_string(data.distance) + "_p" + std::to_string((int)(data.noise * 10000 + 0.5));
        double shot_reps = (double)data.shots.size() * (double)reps;
        const bool probe_match_paths =
            std::getenv("PYMATCHING_PROBE_MATCH_PATHS") != nullptr &&
            std::string(std::getenv("PYMATCHING_PROBE_MATCH_PATHS")) == "1";
        const bool only_probe_match_paths =
            std::getenv("PYMATCHING_ONLY_PROBE_MATCH_PATHS") != nullptr &&
            std::string(std::getenv("PYMATCHING_ONLY_PROBE_MATCH_PATHS")) == "1";
        // Current policy defaults to match-path-selected bounds.
        // Set PYMATCHING_CLUSTER_FROM_MATCH_PATHS=0 only when supplying explicit
        // PYMATCHING_CLUSTER_BOUNDS for a controlled debugging run.
        const bool cluster_from_match_paths = !(
            std::getenv("PYMATCHING_CLUSTER_FROM_MATCH_PATHS") != nullptr &&
            std::string(std::getenv("PYMATCHING_CLUSTER_FROM_MATCH_PATHS")) == "0");
        const bool parameter_schedule_requested =
            env_flag_enabled("PYMATCHING_CLUSTER_PATH_PARAMETER_SCHEDULE");
        const bool cluster_only_mode = env_flag_enabled("PYMATCHING_CLUSTER_ONLY");
        MatchingPathDistanceStats matching_path_stats;
        bool have_matching_path_stats = false;
        size_t global_mistakes = 0;
        std::vector<uint64_t> global_observable_predictions;
        if (!cluster_only_mode) {
            global_observable_predictions = compute_global_observable_predictions(data.graph, data.shots);
            for (size_t k = 0; k < data.shots.size() && k < global_observable_predictions.size(); k++) {
                if (global_observable_predictions[k] != data.shots[k].obs_mask) {
                    global_mistakes++;
                }
            }
            Totals g;
            measure_global_event_counts_only(data.graph, data.shots, reps, g);
            g.mistakes = global_mistakes;
            print_parallel_row(case_name, "global", g, shot_reps);
        }
        if (probe_match_paths) {
            size_t max_probe_shots = 0;
            if (const char* e = std::getenv("PYMATCHING_MATCH_PATH_MAX_SHOTS")) {
                max_probe_shots = static_cast<size_t>(std::stoull(e));
            }
            matching_path_stats = measure_matching_path_distances(data.graph, data.shots, max_probe_shots);
            have_matching_path_stats = true;
            print_matching_path_stats(case_name, matching_path_stats);
        }
        if (only_probe_match_paths) {
            continue;
        }
        size_t actual_workers = 1;
        if (const char* workers_env = std::getenv("PYMATCHING_PARALLEL_WORKERS")) {
            actual_workers = std::max<size_t>(1, static_cast<size_t>(std::stoull(workers_env)));
        }
        // Legacy beta/gamma/lambda automatic hierarchy generation is not used by
        // the current benchmark policy.  The values remain only to initialize the
        // config object when explicit bounds are supplied for debugging.
        const double beta = 2.0e7;
        const double gamma = 2.5e7;
        const double lambda = 2.0;
        std::vector<pm::ProcessingClusterBounds> explicit_bounds;
        if (const char* e = std::getenv("PYMATCHING_CLUSTER_BOUNDS")) {
            // Explicit bounds are kept as a narrow debugging escape hatch.
            explicit_bounds = parse_explicit_cluster_bounds(e);
        } else if (cluster_from_match_paths) {
            const bool parameter_schedule = parameter_schedule_requested;
            pm::cumulative_time_int unit = 0;
            if (parameter_schedule) {
                unit = detector_graph_max_edge_weight_plus_one(data.graph);
            } else if (const char* origin_env = std::getenv("PYMATCHING_CLUSTER_PATH_L_ORIGIN")) {
                unit = static_cast<pm::cumulative_time_int>(std::stoll(origin_env));
            } else if (const char* unit_env = std::getenv("PYMATCHING_CLUSTER_PATH_UNIT")) {
                // Backwards-compatible alias.  Under the current policy this is the
                // unscaled L_origin, not the already-scaled first diameter bound.
                unit = static_cast<pm::cumulative_time_int>(std::stoll(unit_env));
            } else {
                if (!have_matching_path_stats) {
                    size_t max_probe_shots = 0;
                    if (const char* e = std::getenv("PYMATCHING_MATCH_PATH_MAX_SHOTS")) {
                        max_probe_shots = static_cast<size_t>(std::stoull(e));
                    }
                    matching_path_stats = measure_matching_path_distances(data.graph, path_stats_shots, max_probe_shots);
                    have_matching_path_stats = true;
                    print_matching_path_stats(case_name, matching_path_stats);
                }
                unit = choose_path_distance_unit(matching_path_stats, nullptr);
            }
            if (!(unit > 0)) {
                throw std::runtime_error("PYMATCHING_CLUSTER_PATH_L_ORIGIN/PYMATCHING_CLUSTER_PATH_UNIT must be positive.");
            }
            size_t candidate_count = 512;
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_CANDIDATES")) {
                candidate_count = static_cast<size_t>(std::stoull(e));
            }
            size_t max_selected_levels = 64;
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_LEVELS")) {
                max_selected_levels = static_cast<size_t>(std::stoull(e));
            }
            double min_buffer_ratio = 2.2;
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_MIN_BUFFER_RATIO")) {
                min_buffer_ratio = std::stod(e);
            }
            long double path_start_ratio = 1.0L;
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_L_RATIO")) {
                path_start_ratio = std::stold(e);
            }
            if (!(path_start_ratio > 0) || !std::isfinite(static_cast<double>(path_start_ratio))) {
                throw std::runtime_error("PYMATCHING_CLUSTER_PATH_L_RATIO must be positive and finite.");
            }
            bool enforce_legacy_gap_filter =
                env_flag_enabled("PYMATCHING_CLUSTER_PATH_LEGACY_GAP_FILTER");
            bool min_buffer_ratio_first_level_only =
                env_flag_enabled("PYMATCHING_CLUSTER_PATH_MIN_BUFFER_RATIO_FIRST_LEVEL_ONLY");
            bool direct_min_d =
                env_flag_enabled("PYMATCHING_CLUSTER_PATH_DIRECT_MIN_D");
            bool optimize_final_d =
                env_flag_enabled("PYMATCHING_CLUSTER_PATH_OPTIMIZE_FINAL_D");
            bool phi_sequence =
                env_flag_enabled("PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE");
            int path_policy_count = (direct_min_d ? 1 : 0) + (optimize_final_d ? 1 : 0) +
                (phi_sequence ? 1 : 0) + (parameter_schedule ? 1 : 0);
            if (path_policy_count > 1) {
                throw std::runtime_error(
                    "Set at most one of PYMATCHING_CLUSTER_PATH_DIRECT_MIN_D=1, "
                    "PYMATCHING_CLUSTER_PATH_OPTIMIZE_FINAL_D=1, and "
                    "PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE=1, and "
                    "PYMATCHING_CLUSTER_PATH_PARAMETER_SCHEDULE=1.");
            }
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_OPTIMIZE_TARGET_LEVEL")) {
                max_selected_levels = static_cast<size_t>(std::stoull(e));
            }
            bool phi_floor_budget =
                env_flag_enabled("PYMATCHING_CLUSTER_PATH_PHI_FLOOR_BUDGET");
            double phi_floor = 0.5;
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_PHI_FLOOR")) {
                phi_floor = std::stod(e);
            }
            double phi_budget_per_level = 0.24;
            bool have_explicit_phi_budget = false;
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_PHI_BUDGET_PER_LEVEL")) {
                phi_budget_per_level = std::stod(e);
                have_explicit_phi_budget = true;
            }
            if (phi_floor_budget && !have_explicit_phi_budget) {
                if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_PHI_BUDGET_LEVELS")) {
                    size_t budget_levels = static_cast<size_t>(std::stoull(e));
                    if (budget_levels > 1) {
                        phi_budget_per_level =
                            0.98 * (1.0 - phi_floor) / static_cast<double>(budget_levels - 1);
                    }
                } else {
                    phi_budget_per_level = std::min(0.24, 0.48 * (1.0 - phi_floor));
                }
            }
            double optimize_phi_budget_total_fraction = 0.98;
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_OPTIMIZE_PHI_BUDGET_TOTAL_FRACTION")) {
                optimize_phi_budget_total_fraction = std::stod(e);
            }
            size_t optimize_search_rounds = 8;
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_OPTIMIZE_SEARCH_ROUNDS")) {
                optimize_search_rounds = static_cast<size_t>(std::stoull(e));
            }

            double phi_sequence_floor = phi_floor;
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE_FLOOR")) {
                phi_sequence_floor = std::stod(e);
            }
            double phi_sequence_phi2 = 0.0625;
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE_PHI2")) {
                phi_sequence_phi2 = std::stod(e);
            }
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE_B1_RATIO")) {
                double rho = std::stod(e);
                if (!(rho > 1.0) || !std::isfinite(rho)) {
                    throw std::runtime_error(
                        "PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE_B1_RATIO must be finite and greater than 1.");
                }
                phi_sequence_phi2 = (rho - 1.0) / (rho + 2.0);
            }
            double phi_sequence_q = 0.5;
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_PHI_SEQUENCE_Q")) {
                phi_sequence_q = std::stod(e);
            }
            if (optimize_final_d && max_selected_levels == 0) {
                throw std::runtime_error(
                    "PYMATCHING_CLUSTER_PATH_OPTIMIZE_FINAL_D requires a positive target level. "
                    "Set PYMATCHING_CLUSTER_PATH_LEVELS or PYMATCHING_CLUSTER_PATH_OPTIMIZE_TARGET_LEVEL.");
            }
            if (parameter_schedule) {
                explicit_bounds = pm::make_processing_cluster_parameter_schedule_bounds(
                    unit, max_selected_levels, phi_sequence_floor, phi_sequence_q);
            } else if (phi_sequence) {
                auto first_d = ceil_positive_scaled_distance(unit, path_start_ratio);
                explicit_bounds = pm::make_processing_cluster_disjoint_stopping_lemma_bounds_phi_sequence(
                    first_d, max_selected_levels, 1, min_buffer_ratio,
                    min_buffer_ratio_first_level_only, phi_sequence_phi2, phi_sequence_floor,
                    phi_sequence_q);
            } else if (direct_min_d) {
                auto first_d = ceil_positive_scaled_distance(unit, path_start_ratio);
                explicit_bounds = pm::make_processing_cluster_disjoint_stopping_lemma_bounds_direct_min_d(
                    first_d, max_selected_levels, 1, min_buffer_ratio,
                    min_buffer_ratio_first_level_only, phi_floor_budget, phi_floor, phi_budget_per_level);
            } else if (optimize_final_d) {
                auto first_d = ceil_positive_scaled_distance(unit, path_start_ratio);
                explicit_bounds = pm::make_processing_cluster_disjoint_stopping_lemma_bounds_optimize_final_d(
                    first_d, max_selected_levels, 1, min_buffer_ratio,
                    min_buffer_ratio_first_level_only, phi_floor, optimize_phi_budget_total_fraction,
                    optimize_search_rounds);
            } else {
                auto candidates = make_linear_diameter_candidates(unit, path_start_ratio, candidate_count);
                explicit_bounds = pm::make_processing_cluster_disjoint_stopping_lemma_bounds_from_candidates(
                    candidates, max_selected_levels, 1, min_buffer_ratio, enforce_legacy_gap_filter,
                    min_buffer_ratio_first_level_only, phi_floor_budget, phi_floor, phi_budget_per_level);
            }
            if (explicit_bounds.empty()) {
                throw std::runtime_error(
                    parameter_schedule
                        ? "The Algorithm-style parameter schedule produced no usable cluster bounds."
                        : (phi_sequence
                        ? "The phi-sequence matching-path policy produced no usable cluster bounds."
                        : (direct_min_d
                            ? "The direct-min-d matching-path policy produced no usable cluster bounds."
                            : (optimize_final_d
                                ? "The target-level final-d optimizer produced no usable cluster bounds."
                                : "The measured matching-path candidate list produced no usable cluster bounds."))));
            }
            std::cerr << "PATH_SELECTED_BOUNDS_LEGACY_GAP_FILTER "
                      << (enforce_legacy_gap_filter ? 1 : 0) << "\n";
            std::cerr << "PATH_SELECTED_BOUNDS_MIN_BUFFER_RATIO_FIRST_LEVEL_ONLY "
                      << (min_buffer_ratio_first_level_only ? 1 : 0) << "\n";
            std::cerr << "PATH_SELECTED_BOUNDS_DIRECT_MIN_D "
                      << (direct_min_d ? 1 : 0) << "\n";
            std::cerr << "PATH_SELECTED_BOUNDS_OPTIMIZE_FINAL_D "
                      << (optimize_final_d ? 1 : 0) << "\n";
            std::cerr << "PATH_SELECTED_BOUNDS_PHI_SEQUENCE "
                      << (phi_sequence ? 1 : 0) << "\n";
            std::cerr << "PATH_SELECTED_BOUNDS_PARAMETER_SCHEDULE "
                      << (parameter_schedule ? 1 : 0) << "\n";
            if (parameter_schedule) {
                std::cerr << "PATH_SELECTED_BOUNDS_GRAPH_W_MAX " << (unit - 1) << "\n";
                std::cerr << "PATH_SELECTED_BOUNDS_PARAMETER_D1 " << unit << "\n";
                std::cerr << "PATH_SELECTED_BOUNDS_PARAMETER_PHI_FLOOR "
                          << phi_sequence_floor << "\n";
                std::cerr << "PATH_SELECTED_BOUNDS_PARAMETER_Q "
                          << phi_sequence_q << "\n";
            }
            if (phi_sequence) {
                std::cerr << "PATH_SELECTED_BOUNDS_PHI_SEQUENCE_PHI2 "
                          << phi_sequence_phi2 << "\n";
                std::cerr << "PATH_SELECTED_BOUNDS_PHI_SEQUENCE_FLOOR "
                          << phi_sequence_floor << "\n";
                std::cerr << "PATH_SELECTED_BOUNDS_PHI_SEQUENCE_Q "
                          << phi_sequence_q << "\n";
            }
            if (optimize_final_d) {
                std::cerr << "PATH_SELECTED_BOUNDS_OPTIMIZE_TARGET_LEVEL "
                          << max_selected_levels << "\n";
                std::cerr << "PATH_SELECTED_BOUNDS_OPTIMIZE_PHI_BUDGET_TOTAL_FRACTION "
                          << optimize_phi_budget_total_fraction << "\n";
                std::cerr << "PATH_SELECTED_BOUNDS_OPTIMIZE_SEARCH_ROUNDS "
                          << optimize_search_rounds << "\n";
            }
            std::cerr << "PATH_SELECTED_BOUNDS_PHI_FLOOR_BUDGET "
                      << (phi_floor_budget ? 1 : 0) << "\n";
            if (phi_floor_budget) {
                std::cerr << "PATH_SELECTED_BOUNDS_PHI_FLOOR " << phi_floor << "\n";
                std::cerr << "PATH_SELECTED_BOUNDS_PHI_BUDGET_PER_LEVEL "
                          << phi_budget_per_level << "\n";
            }
            print_selected_bounds(
                case_name, unit, parameter_schedule ? 1.0L : path_start_ratio,
                (direct_min_d || optimize_final_d || phi_sequence || parameter_schedule) ? 0 : candidate_count,
                parameter_schedule ? 1.0 : min_buffer_ratio,
                parameter_schedule ? false : min_buffer_ratio_first_level_only,
                direct_min_d, optimize_final_d, phi_sequence, parameter_schedule,
                phi_floor_budget, phi_floor,
                phi_budget_per_level, explicit_bounds);
        } else {
            throw std::runtime_error(
                "Current policy requires match-path selected bounds. Set "
                "PYMATCHING_CLUSTER_FROM_MATCH_PATHS=1 with PYMATCHING_CLUSTER_PATH_L_ORIGIN "
                "or PYMATCHING_CLUSTER_PATH_UNIT, or provide PYMATCHING_CLUSTER_BOUNDS "
                "for debugging.");
        }
        const auto* explicit_bounds_ptr = explicit_bounds.empty() ? nullptr : &explicit_bounds;
        const std::vector<size_t>* max_active_detector_caps_ptr = nullptr;
        if (env_flag_enabled("PYMATCHING_GROWING_ONLY_CLUSTERING")) {
            pm::ProcessingClusterConfig growing_config;
            growing_config.beta = beta;
            growing_config.gamma = gamma;
            growing_config.lambda = lambda;
            growing_config.max_level = explicit_bounds.empty() ? 64 : explicit_bounds.size();
            if (!explicit_bounds.empty()) {
                growing_config.explicit_bounds_by_level = explicit_bounds;
                growing_config.enforce_stopping_lemma_bounds =
                    env_flag_enabled("PYMATCHING_ENFORCE_STOPPING_LEMMA_BOUNDS");
            }
            apply_adaptive_cluster_bounds_env(growing_config);
            size_t max_growing_levels = explicit_bounds.empty() ? 1 : explicit_bounds.size();
            if (const char* e = std::getenv("PYMATCHING_GROWING_ONLY_LEVELS")) {
                max_growing_levels = static_cast<size_t>(std::stoull(e));
            }
            run_and_print_growing_only_benchmark(
                case_name, data.graph, data.shots, reps, growing_config, max_growing_levels);
            if (env_flag_enabled("PYMATCHING_ONLY_GROWING_ONLY_CLUSTERING")) {
                continue;
            }
        }
        if (cluster_only_mode) {
            pm::ProcessingClusterConfig cluster_only_config;
            cluster_only_config.beta = beta;
            cluster_only_config.gamma = gamma;
            cluster_only_config.lambda = lambda;
            cluster_only_config.max_level = explicit_bounds.empty() ? 64 : explicit_bounds.size();
            if (!explicit_bounds.empty()) {
                cluster_only_config.explicit_bounds_by_level = explicit_bounds;
                cluster_only_config.enforce_stopping_lemma_bounds =
                    env_flag_enabled("PYMATCHING_ENFORCE_STOPPING_LEMMA_BOUNDS");
            }
            apply_adaptive_cluster_bounds_env(cluster_only_config);
            auto cluster_counts = run_parallel_cluster_only(data.graph, data.shots, reps, cluster_only_config);
            std::ostringstream mode;
            mode << "parallel_cluster_only_serial";
            if (explicit_bounds_ptr != nullptr) {
                mode << "_explicit_bounds" << explicit_bounds.size();
            }
            if (parameter_schedule_requested) {
                mode << "_parameter_schedule_wmax";
            } else if (cluster_from_match_paths) {
                mode << "_match_path_selected";
            }
            mode << "_lratio";
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_L_RATIO")) {
                mode << e;
            } else {
                mode << "1";
            }
            if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_MIN_BUFFER_RATIO")) {
                mode << "_minbuf" << e;
            }
            print_parallel_row(case_name, mode.str(), cluster_counts, shot_reps);
            continue;
        }
        if (env_flag_enabled("PYMATCHING_EVENT_COUNT_ONLY")) {
            // Default event-count reporting uses the real scheduler-derived
            // imported/rescheduled worker states. The reported ideal event count is
            // now level-by-level critical-path accounting: for each upper level,
            // notifications at or before the previous level's maximum algorithmic
            // stop time are treated as overlappable and not added. The old fresh-run
            // path is kept only behind an explicit debug escape hatch.
            if (env_flag_enabled("PYMATCHING_ALLOW_FRESH_EVENT_COUNT_ONLY")) {
                pm::ProcessingClusterConfig event_count_config;
                event_count_config.beta = beta;
                event_count_config.gamma = gamma;
                event_count_config.lambda = lambda;
                event_count_config.max_level = explicit_bounds.empty() ? 64 : explicit_bounds.size();
                if (!explicit_bounds.empty()) {
                    event_count_config.explicit_bounds_by_level = explicit_bounds;
                    event_count_config.enforce_stopping_lemma_bounds =
                        env_flag_enabled("PYMATCHING_ENFORCE_STOPPING_LEMMA_BOUNDS");
                }
                apply_adaptive_cluster_bounds_env(event_count_config);
                auto event_counts = run_parallel_event_count_only(data.graph, data.shots, reps, event_count_config);
                std::ostringstream mode;
                mode << "parallel_fresh_event_count_only_debug";
                if (explicit_bounds_ptr != nullptr) {
                    mode << "_explicit_bounds" << explicit_bounds.size();
                }
                if (parameter_schedule_requested) {
                    mode << "_parameter_schedule_wmax";
                } else if (cluster_from_match_paths) {
                    mode << "_match_path_selected";
                }
                mode << "_lratio";
                if (const char* e = std::getenv("PYMATCHING_CLUSTER_PATH_L_RATIO")) {
                    mode << e;
                } else {
                    mode << "1";
                }
                print_parallel_row(case_name, mode.str(), event_counts, shot_reps);
                continue;
            }
            auto event_counts = run_parallel(
                data.graph,
                data.shots,
                reps,
                actual_workers,
                64,
                beta,
                gamma,
                lambda,
                true,
                &global_observable_predictions,
                explicit_bounds_ptr,
                max_active_detector_caps_ptr,
                &case_name,
                &data);
            std::ostringstream mode;
            mode << "parallel_scheduler_level_critical_event_counts";
            mode << "_workers" << actual_workers;
            if (explicit_bounds_ptr != nullptr) {
                mode << "_explicit_bounds" << explicit_bounds.size();
            }
            if (parameter_schedule_requested) {
                mode << "_parameter_schedule_wmax";
            } else if (cluster_from_match_paths) {
                mode << "_match_path_selected";
            }
            mode << "_parentless_broadcast_import_growing_only_clusters_level_pipeline";
            print_parallel_row(case_name, mode.str(), event_counts, shot_reps);
            continue;
        }
        auto p1 = run_parallel(
            data.graph,
            data.shots,
            reps,
            actual_workers,
            64,
            beta,
            gamma,
            lambda,
            true,
            &global_observable_predictions,
            explicit_bounds_ptr,
            max_active_detector_caps_ptr,
            &case_name,
            &data);
        std::ostringstream mode;
        mode << "parallel_lockstep_full_graph_workers";
        mode << "_workers" << actual_workers;
        if (explicit_bounds_ptr != nullptr) {
            mode << "_explicit_bounds" << explicit_bounds.size();
        } else {
            mode << "_idealw64_b" << beta << "_g" << gamma << "_l" << lambda;
        }
        if (parameter_schedule_requested) {
            mode << "_parameter_schedule_wmax";
        } else if (cluster_from_match_paths) {
            mode << "_match_path_selected";
        }
        mode << "_parentless_broadcast_import_growing_only_clusters_level_pipeline";
        print_parallel_row(case_name, mode.str(), p1, shot_reps);
    }
}
