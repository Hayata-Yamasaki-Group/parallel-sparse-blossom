#ifndef STIM_STUB_H
#define STIM_STUB_H
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
namespace stim {
struct DemTarget {
    bool is_relative_detector_id() const { return false; }
    bool is_observable_id() const { return false; }
    bool is_separator() const { return false; }
    size_t val() const { return 0; }
    size_t raw_id() const { return 0; }
};
struct DemInstruction {
    std::vector<double> arg_data;
    std::vector<DemTarget> target_data;
};
class DetectorErrorModel {
public:
    size_t count_detectors() const { return 0; }
    size_t count_observables() const { return 0; }
    template <typename F>
    void iter_flatten_error_instructions(const F&) const {}
};
}
#endif
