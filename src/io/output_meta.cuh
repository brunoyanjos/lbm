#pragma once
#include <string>
#include "../app/cuda_config.cuh"
#include "../app/benchmark.cuh"

namespace io
{
    void write_performance(const std::string &out_dir,
                           const CudaConfig &cfg,
                           const app::BenchmarkResult &r);
}
