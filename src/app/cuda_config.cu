#include "cuda_config.cuh"
#include "../lbm/lbm_tuner.cuh"
#include "../core/geometry.h"
#include "../lbm/stencil_active.cuh"
#include <iostream>

CudaConfig make_config(int nx, int ny)
{
    CudaConfig cfg;

    cfg.block = dim3(16, 16);

    cfg.grid = dim3(
        (nx + cfg.block.x - 1) / cfg.block.x,
        (ny + cfg.block.y - 1) / cfg.block.y);

    cfg.shared_bytes =
        cfg.block.x * cfg.block.y * (Stencil::Q - 1) * sizeof(real_t);

    return cfg;
}
