#pragma once

#include "../../../core/types.cuh"
#include "../../../geometries/active_geometry.cuh"
#include "../../../core/math_utils.cuh"
#include "../../../core/linear_solver.cuh"
#include "../../stencil_active.cuh"
#include "../../domain/mask_utils.cuh"
#include "../../state/lbm_state.cuh"
#include <cstdint>
#include <cstdio>

#include "system_variables.cuh"
#include "helpers.cuh"

__device__ __forceinline__ void
apply_boundary_dirichlet(
    real_t *__restrict__ pop,
    uint32_t valid_mask,
    real_t &rho,
    real_t ux, real_t uy,
    real_t &mxx, real_t &mxy, real_t &myy)
{
    const uint32_t outgoing_mask = valid_mask;
    const uint32_t incoming_mask = mask_opp(valid_mask);

    DirichletVariable2D S{};

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        real_t cx, cy, Hxx, Hxy, Hyy;
        Stencil::basis2(i, cx, cy, Hxx, Hxy, Hyy);

        const real_t w = Stencil::w(i);
        const real_t Hx_factor = w * Stencil::as2 * cx;
        const real_t Hy_factor = w * Stencil::as2 * cy;
        const real_t Hxx_factor = w * Stencil::as4 * r::half * Hxx;
        const real_t Hxy_factor = w * Stencil::as4 * Hxy;
        const real_t Hyy_factor = w * Stencil::as4 * r::half * Hyy;

        if (dir_valid(incoming_mask, i))
            accumulate_incoming_dirichlet(
                S, pop[i],
                w, Hx_factor, Hy_factor,
                Hxx_factor, Hxy_factor, Hyy_factor,
                Hxx, Hxy, Hyy);

        if (dir_valid(outgoing_mask, i))
            accumulate_outgoing_dirichlet(
                S,
                w, Hx_factor, Hy_factor,
                Hxx_factor, Hxy_factor, Hyy_factor);
    }

    S.normalize_moments();

    real_t Aflat[9]{};
    real_t b[3]{};
    real_t m[3]{};

    build_dirichlet_system(S, ux, uy, Aflat, b);
    gaussianElimination<3>(Aflat, b, m);

    mxx = m[0];
    mxy = m[1];
    myy = m[2];

    rho = recover_dirichlet_rho(S, ux, uy, mxx, mxy, myy);
}