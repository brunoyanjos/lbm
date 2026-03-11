#pragma once

#include "../../../core/types.cuh"
#include "../../../core/linear_solver.cuh"
#include "../../stencil_active.cuh"
#include "../../domain/mask_utils.cuh"
#include "system_variables.cuh"

__device__ __forceinline__ void accumulate_incoming_fluid(
    FluidSystem2D &S,
    real_t pop_i,
    real_t w,
    real_t Hx_factor,
    real_t Hy_factor,
    real_t Hxx_factor,
    real_t Hxy_factor,
    real_t Hyy_factor,
    real_t cx,
    real_t cy,
    real_t Hxx,
    real_t Hxy,
    real_t Hyy)
{
    S.rho_I += pop_i;

    S.ux_I += pop_i * cx;
    S.uy_I += pop_i * cy;

    S.mxx_I += pop_i * Hxx;
    S.mxy_I += pop_i * Hxy;
    S.myy_I += pop_i * Hyy;

    S.ux_rho += w * cx;
    S.ux_ux += Hx_factor * cx;
    S.ux_uy += Hy_factor * cx;
    S.ux_mxx += Hxx_factor * cx;
    S.ux_mxy += Hxy_factor * cx;
    S.ux_myy += Hyy_factor * cx;

    S.uy_rho += w * cy;
    S.uy_ux += Hx_factor * cy;
    S.uy_uy += Hy_factor * cy;
    S.uy_mxx += Hxx_factor * cy;
    S.uy_mxy += Hxy_factor * cy;
    S.uy_myy += Hyy_factor * cy;

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

__device__ __forceinline__ void accumulate_outgoing_fluid(
    FluidSystem2D &S,
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

__device__ __forceinline__ real_t recover_fluid_rho(
    const FluidSystem2D &S,
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