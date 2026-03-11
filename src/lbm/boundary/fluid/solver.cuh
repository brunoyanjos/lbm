#pragma once

#include "../../../core/types.cuh"
#include "../../../core/linear_solver.cuh"
#include "../../domain/mask_utils.cuh"
#include "../../../geometries/active_geometry.cuh"

#include "system_variables.cuh"
#include "helpers.cuh"
#include "non_linear_system.cuh"

#include <cstdint>

__device__ inline void evaluate_fluid_node(
    real_t *__restrict__ pop,
    uint32_t valid_mask,
    real_t &rho,
    real_t &ux, real_t &uy,
    real_t &mxx, real_t &mxy, real_t &myy)
{
    const uint32_t outgoing_mask = valid_mask;
    const uint32_t incoming_mask = mask_opp(valid_mask);

    FluidSystem2D S{};

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        real_t cx, cy, Hxx, Hxy, Hyy;
        Stencil::basis2(i, cx, cy, Hxx, Hxy, Hyy);

        const real_t w = Stencil::w(i);
        const real_t Hx_factor = w * Stencil::as2 * cx;
        const real_t Hy_factor = w * Stencil::as2 * cy;
        const real_t Hxx_factor = w * Stencil::as4 * real_t(0.5) * Hxx;
        const real_t Hxy_factor = w * Stencil::as4 * Hxy;
        const real_t Hyy_factor = w * Stencil::as4 * real_t(0.5) * Hyy;

        if (dir_valid(incoming_mask, i))
            accumulate_incoming_fluid(
                S, pop[i],
                w, Hx_factor, Hy_factor,
                Hxx_factor, Hxy_factor, Hyy_factor,
                cx, cy,
                Hxx, Hxy, Hyy);

        if (dir_valid(outgoing_mask, i))
            accumulate_outgoing_fluid(
                S,
                w, Hx_factor, Hy_factor,
                Hxx_factor, Hxy_factor, Hyy_factor);
    }

    S.normalize_incoming();

    FluidNewtonCoefficients2D C{};
    build_fluid_newton_coefficients(S, C);

    solve_fluid_newton(C, ux, uy, mxx, mxy, myy);

    rho = recover_fluid_rho(S, ux, uy, mxx, mxy, myy);
}