#pragma once

#include "../../lbm/state/lbm_state.cuh"
#include "../../core/types.cuh"

#include <string>
#include <cstdint>
#include <cstddef>
#include <vector>

namespace io
{
    struct FlowAverages
    {
        std::int64_t samples = 0;
        std::vector<double> ux;
        std::vector<double> uy;
        std::vector<double> uxux;
        std::vector<double> uxuy;
        std::vector<double> uyuy;
    };

    [[nodiscard]] FlowAverages make_flow_averages(std::size_t node_count);

    real_t compute_ke_host_2d(const LBMState &state);

    real_t compute_ke_and_accumulate_flow_averages_host_2d(const LBMState &state,
                                                           FlowAverages &averages);

    void tke_bin_append(const std::string &out_dir, int t, double ke);

    void seed_tke_history_from_checkpoint(const std::string &checkpoint_path_or_dir,
                                          const std::string &out_dir,
                                          int checkpoint_step);

    bool seed_flow_averages_from_checkpoint(const std::string &checkpoint_path_or_dir,
                                            FlowAverages &averages,
                                            int checkpoint_step);

    void write_flow_averages(const FlowAverages &averages,
                             int step,
                             const std::string &out_dir);

    void write_flow_averages_checkpoint(const FlowAverages &averages,
                                        int step,
                                        const std::string &out_dir);
}
