#include "cuda_config.cuh"
#include "../lbm/lbm_tuner.cuh"
#include "../core/active_geometry.cuh"
#include "../lbm/stencil_active.cuh"
#include <iostream>

CudaConfig make_config(
#if defined(LBM_STENCIL_D2V17)
    int shared_memory
#endif
)
{
    CudaConfig cfg;

#if defined(LBM_STENCIL_D2Q9)
    cfg.block = dim3(32, 32);
#elif defined(LBM_STENCIL_D2V17)
    size_t bytes_per_lattice = (Stencil::Q - 1) * sizeof(real_t);
    size_t max_lattices = shared_memory / bytes_per_lattice;

    cfg.block = find_optimal_block(max_lattices);
#else
    cfg.block = dim3(32, 32);
#endif

    cfg.grid = dim3(
        (NX + cfg.block.x - 1) / cfg.block.x,
        (NY + cfg.block.y - 1) / cfg.block.y);

    cfg.shared_bytes =
        cfg.block.x * cfg.block.y * (Stencil::Q - 1) * sizeof(real_t);

    return cfg;
}
