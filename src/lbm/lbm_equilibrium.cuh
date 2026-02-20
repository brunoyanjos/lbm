#pragma once

#include "../core/types.cuh"
#include "stencil_active.cuh"

__device__ __forceinline__ void equilibrium(real_t *pop, real_t rho, real_t ux, real_t uy)
{
#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const real_t cx = static_cast<real_t>(Stencil::cx(i));
        const real_t cy = static_cast<real_t>(Stencil::cy(i));

        const real_t Hxx = cx * cx - Stencil::cs2;
        const real_t Hxy = cx * cy;
        const real_t Hyy = cy * cy - Stencil::cs2;

        pop[i] = Stencil::w(i) * rho *
                 (static_cast<real_t>(1) + Stencil::as2 * ux * cx + Stencil::as2 * uy * cy +
                  (Stencil::as4 * static_cast<real_t>(0.5)) * ux * ux * Hxx +
                  (Stencil::as4 * static_cast<real_t>(0.5)) * uy * uy * Hyy +
                  Stencil::as4 * ux * uy * Hxy);
    }
}
