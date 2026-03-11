#pragma once

#include "../../../core/types.cuh"
#include "../../../core/linear_solver.cuh"
#include "../../stencil_active.cuh"
#include "../../domain/mask_utils.cuh"
#include "system_variables.cuh"

__device__ __forceinline__ void accumulate_incoming_dirichlet(
    DirichletVariable2D &S,
    real_t pop_i,
    real_t w,
    real_t Hx_factor,
    real_t Hy_factor,
    real_t Hxx_factor,
    real_t Hxy_factor,
    real_t Hyy_factor,
    real_t Hxx,
    real_t Hxy,
    real_t Hyy)
{
    S.rho_I += pop_i;

    S.mxx_I += pop_i * Hxx;
    S.mxy_I += pop_i * Hxy;
    S.myy_I += pop_i * Hyy;

    S.mxx_rho += w * Hxx;
    S.mxx_ux += Hx_factor * Hxx;
    S.mxx_uy += Hy_factor * Hxx;
    S.mxx_mxx += Hxx_factor * Hxx;
    S.mxx_mxy += Hxy_factor * Hxx;
    S.mxx_myy += Hyy_factor * Hxx;

    S.mxy_rho += w * Hxy;
    S.mxy_ux += Hx_factor * Hxy;
    S.mxy_uy += Hy_factor * Hxy;
    S.mxy_mxx += Hxx_factor * Hxy;
    S.mxy_mxy += Hxy_factor * Hxy;
    S.mxy_myy += Hyy_factor * Hxy;

    S.myy_rho += w * Hyy;
    S.myy_ux += Hx_factor * Hyy;
    S.myy_uy += Hy_factor * Hyy;
    S.myy_mxx += Hxx_factor * Hyy;
    S.myy_mxy += Hxy_factor * Hyy;
    S.myy_myy += Hyy_factor * Hyy;
}

__device__ __forceinline__ void accumulate_outgoing_dirichlet(
    DirichletVariable2D &S,
    real_t w,
    real_t Hx_factor,
    real_t Hy_factor,
    real_t Hxx_factor,
    real_t Hxy_factor,
    real_t Hyy_factor)
{
    S.rho_rho += w;
    S.rho_ux += Hx_factor;
    S.rho_uy += Hy_factor;
    S.rho_mxx += (r::one - Geometry::OMEGA) * Hxx_factor;
    S.rho_mxy += (r::one - Geometry::OMEGA) * Hxy_factor;
    S.rho_myy += (r::one - Geometry::OMEGA) * Hyy_factor;
    S.rho_uxux += Geometry::OMEGA * Hxx_factor;
    S.rho_uxuy += Geometry::OMEGA * Hxy_factor;
    S.rho_uyuy += Geometry::OMEGA * Hyy_factor;
}

__device__ __forceinline__ void build_dirichlet_system(
    const DirichletVariable2D &S,
    real_t ux,
    real_t uy,
    real_t *Aflat,
    real_t *b)
{
    auto A = [&](int r, int c) -> real_t &
    { return Aflat[r * 3 + c]; };

    const real_t rho_part =
        S.rho_rho +
        ux * S.rho_ux + uy * S.rho_uy +
        ux * ux * S.rho_uxux + ux * uy * S.rho_uxuy + uy * uy * S.rho_uyuy;

    A(0, 0) = S.mxx_mxx - S.rho_mxx * S.mxx_I;
    A(0, 1) = S.mxx_mxy - S.rho_mxy * S.mxx_I;
    A(0, 2) = S.mxx_myy - S.rho_myy * S.mxx_I;
    b[0] = S.mxx_I * rho_part - (S.mxx_rho + ux * S.mxx_ux + uy * S.mxx_uy);

    A(1, 0) = S.mxy_mxx - S.rho_mxx * S.mxy_I;
    A(1, 1) = S.mxy_mxy - S.rho_mxy * S.mxy_I;
    A(1, 2) = S.mxy_myy - S.rho_myy * S.mxy_I;
    b[1] = S.mxy_I * rho_part - (S.mxy_rho + ux * S.mxy_ux + uy * S.mxy_uy);

    A(2, 0) = S.myy_mxx - S.rho_mxx * S.myy_I;
    A(2, 1) = S.myy_mxy - S.rho_mxy * S.myy_I;
    A(2, 2) = S.myy_myy - S.rho_myy * S.myy_I;
    b[2] = S.myy_I * rho_part - (S.myy_rho + ux * S.myy_ux + uy * S.myy_uy);
}

__device__ __forceinline__ real_t recover_dirichlet_rho(
    const DirichletVariable2D &S,
    real_t ux,
    real_t uy,
    real_t mxx,
    real_t mxy,
    real_t myy)
{
    const real_t rho_denominator =
        S.rho_rho + ux * S.rho_ux + uy * S.rho_uy +
        mxx * S.rho_mxx + mxy * S.rho_mxy + myy * S.rho_myy +
        ux * ux * S.rho_uxux + ux * uy * S.rho_uxuy + uy * uy * S.rho_uyuy;

    return S.rho_I / rho_denominator;
}