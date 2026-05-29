#pragma once
#include <cuda_runtime.h>
#include <iostream>
#include <vector>

#include "../core/geometry.h"
#include "../core/local_domain.cuh"
#include "../core/physics.h"
#include "../core/simulation_config.h"
#include "../app/cuda_config.cuh"
#include "../app/domain_partition.cuh"
#include "../lbm/stencil_active.cuh"

inline void print_simulation_summary(const cudaDeviceProp &prop,
                                     const std::vector<app::DomainPartition> &partitions)
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
#elif defined(LBM_STENCIL_D2V37)
    std::cout << "D2V37\n";
#else
    std::cout << "UNKNOWN\n";
#endif

    std::cout << "Q                 : " << Stencil::Q << "\n";
    std::cout << "Stencil radius    : " << Stencil::radius << "\n";
    std::cout << "cs^2              : " << Stencil::cs2 << "\n\n";

    // Regularization
    std::cout << "Regularization    : order " << REG_ORDER << "\n";
    std::cout << "Recurrence        : " << (USE_RECURRENCE ? "on" : "off") << "\n\n";

    // Domain
    std::cout << "Domain size       : "
              << NX << " x " << NY << "\n";
    std::cout << "Total nodes       : "
              << NX * NY << "\n\n";

    // Partitions
    std::cout << "Partitions        : " << partitions.size() << "\n";
    for (size_t i = 0; i < partitions.size(); ++i)
    {
        const app::DomainPartition &partition = partitions[i];
        const LocalDomain local = make_local_domain(partition.y_begin,
                                                    partition.y_end,
                                                    partition.halo);
        std::cout << "Partition " << i << "       : device=" << partition.device_id
                  << " y=[" << partition.y_begin << "," << partition.y_end << ")"
                  << " local_ny=" << partition.local_ny
                  << " halo=" << partition.halo
                  << " storage_ny=" << local.storage_ny
                  << " storage_nodes=" << local.N << "\n";
    }
    std::cout << "\n";

    // Physics / nondimensional
    std::cout << "Re                : " << RE << "\n";
    std::cout << "U_lid             : " << U_LID << "\n";
    std::cout << "L_char            : " << L_CHAR << "\n";
    std::cout << "nu (visc)         : " << VISC << "\n";
    std::cout << "tau               : " << TAU << "\n";
    std::cout << "omega             : " << OMEGA << "\n";
    std::cout << "rho_0             : " << RHO_0 << "\n\n";

    std::cout << "====================================================\n\n";
}
