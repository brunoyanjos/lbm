#pragma once
#include <string>
#include "cuda_config.cuh"
#include "run_context.cuh"

namespace app
{
    void run(const CudaConfig &cfg, const RunContext &ctx);
}