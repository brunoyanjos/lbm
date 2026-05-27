#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
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

    void validate_flow_averages_shape(const FlowAverages &averages);

    void sample_flow_averages_node(FlowAverages &averages,
                                   std::size_t idx,
                                   double ux,
                                   double uy,
                                   double next_sample);

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
