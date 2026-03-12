#pragma once

#include "../../../core/types.cuh"
#include "../common/moment_ids.cuh"
#include "../common/nonlinear_ids.cuh"
#include "system.cuh"

__device__ __forceinline__ real_t recover_fluid_rho(
    const FluidSystem2D &S,
    real_t ux,
    real_t uy,
    real_t mxx,
    real_t mxy,
    real_t myy)
{
    const real_t denom =
        S.rho.lin.template get<MomentId2D::rho>() +
        ux * S.rho.lin.template get<MomentId2D::ux>() +
        uy * S.rho.lin.template get<MomentId2D::uy>() +
        mxx * S.rho.lin.template get<MomentId2D::mxx>() +
        mxy * S.rho.lin.template get<MomentId2D::mxy>() +
        myy * S.rho.lin.template get<MomentId2D::myy>() +
        ux * ux * S.rho.nonlin.template get<NonlinearId2D::uxux>() +
        ux * uy * S.rho.nonlin.template get<NonlinearId2D::uxuy>() +
        uy * uy * S.rho.nonlin.template get<NonlinearId2D::uyuy>();

    return S.incomings.template get<MomentId2D::rho>() / denom;
}