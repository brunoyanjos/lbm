#pragma once

#include "core/indexing.cuh"
#include "core/types.cuh"
#include "lbm/interface/normal.cuh"
#include "lbm/state/lbm_state.cuh"
#include "lbm/stencil_active.cuh"

template <int RegOrder, bool HighOrder>
[[nodiscard]]
__device__ __forceinline__ Normal evaluate_normal(int c, int x, int y,
                                                  const LBMStateFor<RegOrder, HighOrder> &S)
{
    Normal n{};

    for (int i = 0; i < Stencil::Q; ++i)
    {
        const int cx = Stencil::cx(i);
        const int cy = Stencil::cy(i);

        const size_t n_idx = idxGlobalPeriodic(x + cx, y + cy);

        const real_t rhoA = S.d_rhoA[c][n_idx];
        const real_t rhoB = S.d_rhoB[c][n_idx];
        const real_t rho = rhoA + rhoB;
        const real_t phi = rho > r::zero ? (rhoA - rhoB) / rho : r::zero;

        n.x += Stencil::w(i) * phi * cx * Stencil::as2;
        n.y += Stencil::w(i) * phi * cy * Stencil::as2;
    }

    const real_t norm = dsqrt(n.x * n.x + n.y * n.y);
    const real_t inv_norm = norm > r::zero ? r::one / norm : r::zero;

    n.x *= inv_norm;
    n.y *= inv_norm;

    return n;
}
