#pragma once

#include "../../../core/types.cuh"
#include "../../../core/math_utils.cuh"
#include "../../../core/linear_solver.cuh"
#include "../../../core/lbm_features.cuh"
#include "../../../geometries/active_geometry.cuh"

#include "../../stencil_active.cuh"
#include "../../domain/mask_utils.cuh"

#include "../common/factors.cuh"
#include "../common/make_factors.cuh"
#include "system.cuh"
#include "accumulate.cuh"
#include "build_system.cuh"
#include "recover_rho.cuh"

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

    DirichletSystem2D S{};

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const auto B = Stencil::basis(i);
        const real_t w = Stencil::w(i);
        const Factors2D F = make_factors(B, w);

        if (dir_valid(incoming_mask, i))
            accumulate_incoming_dirichlet(S, pop[i], B, F);

        if (dir_valid(outgoing_mask, i))
            accumulate_outgoing_dirichlet(S, F);
    }

    S.normalize_known();

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