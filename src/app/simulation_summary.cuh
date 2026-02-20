#pragma once
#include <cuda_runtime.h>
#include <iostream>

#include "../core/geometry.h"
#include "../app/cuda_config.cuh"
#include "../lbm/stencil_active.cuh"

inline void print_simulation_summary(const CudaConfig &cfg,
                                     const cudaDeviceProp &prop)
{
    std::cout << "\n================ Simulation Summary ================\n";

    // GPU
    std::cout << "GPU              : " << prop.name << "\n";
    std::cout << "Compute Capability: "
              << prop.major << "." << prop.minor << "\n";
    std::cout << "Shared mem / block: "
              << prop.sharedMemPerBlock / 1024 << " KB\n\n";

    // Stencil
    std::cout << "Stencil           : ";
#if defined(LBM_STENCIL_D2Q9)
    std::cout << "D2Q9\n";
#elif defined(LBM_STENCIL_D2V17)
    std::cout << "D2V17\n";
#else
    std::cout << "UNKNOWN\n";
#endif

    std::cout << "Q                 : " << Stencil::Q << "\n";
    std::cout << "cs^2              : " << Stencil::cs2 << "\n\n";

    // Domain
    std::cout << "Domain size       : "
              << NX << " x " << NY << "\n";
    std::cout << "Total nodes       : "
              << NX * NY << "\n\n";

    // Kernel config
    const size_t threads_per_block = cfg.block.x * cfg.block.y;

    std::cout << "Block dim         : ("
              << cfg.block.x << ", "
              << cfg.block.y << ", "
              << cfg.block.z << ")\n";

    std::cout << "Grid dim          : ("
              << cfg.grid.x << ", "
              << cfg.grid.y << ", "
              << cfg.grid.z << ")\n";

    std::cout << "Threads / block   : "
              << threads_per_block << "\n";

    std::cout << "Shared bytes/block: "
              << cfg.shared_bytes << " bytes\n";

    std::cout << "====================================================\n\n";
}
