#include "cuda_config.cuh"
#include "../lbm/lbm_tuner.cuh"
#include "../geometries/active_geometry.cuh"
#include "../lbm/stencil_active.cuh"
#include <iostream>

CudaConfig make_config()
{
    CudaConfig cfg;

    cfg.block = dim3(32, 3);

    cfg.grid = dim3(
        (Geometry::NX + cfg.block.x - 1) / cfg.block.x,
        (Geometry::NY + cfg.block.y - 1) / cfg.block.y);

    cfg.shared_bytes =
        cfg.block.x * cfg.block.y * (Stencil::Q - 1) * sizeof(real_t);

    return cfg;
}
