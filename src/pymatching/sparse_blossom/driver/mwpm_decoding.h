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

#ifndef PYMATCHING2_MWPM_DECODING_H
#define PYMATCHING2_MWPM_DECODING_H

#include <cstddef>
#include <string>

#include "pymatching/sparse_blossom/matcher/mwpm.h"
#include "stim.h"

namespace pm {

class UserGraph;

struct ExtendedMatchingResult {
    std::vector<uint8_t> obs_crossed;
    total_weight_int weight;
    ExtendedMatchingResult();
    explicit ExtendedMatchingResult(size_t num_observables);

    bool operator==(const ExtendedMatchingResult& rhs) const;

    bool operator!=(const ExtendedMatchingResult& rhs) const;

    ExtendedMatchingResult(std::vector<uint8_t> obs_crossed, total_weight_int weight);

    void reset();

    ExtendedMatchingResult& operator+=(const ExtendedMatchingResult& rhs);
    ExtendedMatchingResult operator+(const ExtendedMatchingResult& rhs) const;
};

inline void ExtendedMatchingResult::reset() {
    std::fill(obs_crossed.begin(), obs_crossed.end(), 0);
    weight = 0;
}

void fill_bit_vector_from_obs_mask(pm::obs_int obs_mask, uint8_t* obs_begin_ptr, size_t num_observables);
obs_int bit_vector_to_obs_mask(const std::vector<uint8_t>& bit_vector);

enum class DecoderMode : uint8_t {
    SPARSE_BLOSSOM,
    PARALLEL_SPARSE_BLOSSOM,
};

enum class ParallelDecoderRuntimeMode : uint8_t {
    PRODUCTION_SAFE,
    PAPER_SEMANTICS_TEST,
};

struct ParallelSparseBlossomConfig {
    double beta = 72;
    double gamma = 291;
    double lambda = 320;
    size_t max_level = 64;
    size_t num_workers = 1;
    ParallelDecoderRuntimeMode runtime_mode = ParallelDecoderRuntimeMode::PRODUCTION_SAFE;
};

struct DecoderConfig {
    DecoderMode mode = DecoderMode::SPARSE_BLOSSOM;
    bool edge_correlations = false;
    ParallelSparseBlossomConfig parallel_config;
};

Mwpm detector_error_model_to_mwpm(
    const stim::DetectorErrorModel& detector_error_model,
    pm::weight_int num_distinct_weights,
    bool ensure_search_flooder_included = false,
    bool enable_correlations = false);

MatchingResult decode_detection_events_for_up_to_64_observables(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, bool edge_correlations);

/// Initializes an MWPM instance with the supplied detection events without advancing the flooder
/// timeline. This is the setup stage used by correctness-first experimental schedulers.
void initialize_mwpm_for_detection_events(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

/// Advances the internal sparse blossom state until no further MWPM notifications remain.
/// This is the shatter-before-extract stopping point used by experimental hierarchical schedulers.
/// Throws if the remaining state has unmatched alternating-tree structure and therefore no perfect
/// matching exists for the supplied detection events.
void process_timeline_until_completion(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

/// Extracts the compact final matching result from a quiescent MWPM state without replaying
/// the flooder timeline.  This matches the fast global sparse-blossom decode path used when
/// the graph has at most 64 observables: the result is kept as an obs_mask + weight instead
/// of being expanded into an ExtendedMatchingResult vector.  The caller must ensure the MWPM
/// is already quiescent without shattering.
MatchingResult extract_compact_result_from_current_mwpm_state(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

/// Extracts the final matching result from a quiescent MWPM state without replaying the flooder
/// timeline. The caller must ensure the MWPM is already quiescent without shattering.
ExtendedMatchingResult extract_result_from_current_mwpm_state(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

ExtendedMatchingResult decode_detection_events_in_user_graph(
    pm::UserGraph& graph,
    const std::vector<uint64_t>& detection_events,
    const DecoderConfig& config = DecoderConfig());

/// Used to decode detection events for an existing Mwpm object `mwpm', and a vector of
/// detection event indices `detection_events'. The predicted observables are XOR-ed into an
/// existing uint8_t array with at least `mwpm.flooder.graph.num_observables' elements,
/// the pointer to the first element of which is passed as the `obs_begin_ptr' argument.
/// The weight of the MWPM solution is added to the `weight' argument.
void decode_detection_events(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    uint8_t* obs_begin_ptr,
    pm::total_weight_int& weight,
    bool edge_correlations);

/// Decode detection events using a Mwpm object and vector of detection event indices
/// Returns the compressed edges in the matching: the pairs of detection events that are
/// matched to each other via paths.
void decode_detection_events_to_match_edges(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

/// Decode detection events using a Mwpm object and vector of detection event indices.
/// Returns the edges in the matching: these are pairs of *detectors* forming *edges* in the
/// matching solution (rather than pairs of detection *events* matched via *paths* as returned
/// instead by `decode_detection_events_to_match_edges`).
void decode_detection_events_to_edges(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<int64_t>& edges);

void decode_detection_events_to_edges_with_edge_correlations(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<int64_t>& edges);

}  // namespace pm

#endif  // PYMATCHING2_MWPM_DECODING_H
