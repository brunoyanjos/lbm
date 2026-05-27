#pragma once

#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/system_data.cuh"

#include "lbm/boundary/fluid/build_newton_step.cuh"
#include "lbm/boundary/fluid/max_rel_step.cuh"

template <std::size_t N, std::size_t E>
__device__ __forceinline__ void solve_fluid_newton(const SystemData<N, E> &S, NodeMoments &M)
{
    SystemData<N> G{};

    real_t error = r::one;
    int it = 0;
    constexpr int it_max = 50;
    constexpr real_t tol = real_t(1e-6);

    while (error > tol && it++ < it_max)
    {
        const NodeMoments old = M;

        build_newton_step(G, S, M);
        gaussianElimination(G);

        M.ux += G.x[0];
        M.uy += G.x[1];
        M.mxx += G.x[2];
        M.mxy += G.x[3];
        M.myy += G.x[4];

        error = max_rel_step(M, old);
    }
}