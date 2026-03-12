#pragma once

#include "../../../core/types.cuh"
#include "../../../core/linear_solver.cuh"
#include "../../../geometries/active_geometry.cuh"
#include "../../domain/mask_utils.cuh"
#include "../../stencil_active.cuh"

#include "../common/factors.cuh"
#include "../common/make_factors.cuh"

#include "system.cuh"
#include "accumulate.cuh"
#include "newton_coefficients.cuh"
#include "build_newton_system.cuh"
#include "solve_newton.cuh"
#include "recover_rho.cuh"

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
        const auto B = Stencil::basis(i);
        const real_t w = Stencil::w(i);
        const Factors2D F = make_factors(B, w);

        if (dir_valid(incoming_mask, i))
            accumulate_incoming_fluid(S, pop[i], B, F);

        if (dir_valid(outgoing_mask, i))
            accumulate_outgoing_fluid(S, F);
    }

    S.normalize_known();

    FluidNewtonCoefficients2D C{};
    build_fluid_newton_coefficients(S, C);

    solve_fluid_newton(C, ux, uy, mxx, mxy, myy);

    rho = recover_fluid_rho(S, ux, uy, mxx, mxy, myy);
}