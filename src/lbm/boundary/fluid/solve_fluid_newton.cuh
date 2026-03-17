#pragma once

#include "lbm/boundary/common/system_data.cuh"
#include "lbm/moment/node_moments.cuh"

template <std::size_t N, std::size_t E>
__device__ __forceinline__ void solve_fluid_newton(const SystemData<N, E> &S, NodeMoments &u)
{
    SystemData<N> G{};

    real_t error = r::one;
    int it = 0;
    constexpr int it_max = 50;
    constexpr real_t tol = real_t(1e-6);

    while (error > tol && it++ < it_max)
    {
        const NodeMoments old = u;

        build_newton_step(G, S, u);
        gaussianElimination(G);

        u.ux += G.x[0];
        u.uy += G.x[1];
        u.mxx += G.x[2];
        u.mxy += G.x[3];
        u.myy += G.x[4];

        error = max_rel_step(u, old);
    }
}