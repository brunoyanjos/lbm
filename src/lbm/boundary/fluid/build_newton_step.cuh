#pragma once

#include <cstddef>

#include "core/types.cuh"
#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/system_data.cuh"

template <std::size_t N, std::size_t E>
__device__ __forceinline__ real_t eval_row(
    const SystemData<N, E> &S, const NodeMoments &M, int i)
{
    const real_t ux2 = M.ux * M.ux;
    const real_t uxuy = M.ux * M.uy;
    const real_t uy2 = M.uy * M.uy;

    return M.ux * S.coeff(i, 0) + M.uy * S.coeff(i, 1) +
           ux2 * S.coeff(i, 2) + uxuy * S.coeff(i, 3) + uy2 * S.coeff(i, 4) +
           M.mxx * S.coeff(i, 5) + M.mxy * S.coeff(i, 6) + M.myy * S.coeff(i, 7) -
           S.b[i];
}

template <std::size_t N, std::size_t E>
__device__ __forceinline__ void build_newton_step(SystemData<N> &G,
                                                  const SystemData<N, E> &S,
                                                  const NodeMoments &M)
{
#pragma unroll
    for (int i = 0; i < int(N); ++i)
    {
        G.coeff(i, 0) = S.coeff(i, 0) + r::two * S.coeff(i, 2) * M.ux + S.coeff(i, 3) * M.uy;
        G.coeff(i, 1) = S.coeff(i, 1) + S.coeff(i, 3) * M.ux + r::two * S.coeff(i, 4) * M.uy;
        G.coeff(i, 2) = S.coeff(i, 5);
        G.coeff(i, 3) = S.coeff(i, 6);
        G.coeff(i, 4) = S.coeff(i, 7);

        G.b[i] = -eval_row(S, M, i);
    }
}