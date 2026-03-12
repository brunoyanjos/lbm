#pragma once

#include "../../../core/types.cuh"
#include "../common/moment_ids.cuh"
#include "../common/nonlinear_ids.cuh"
#include "system.cuh"

__device__ __forceinline__ void build_fluid_system(
    const FluidSystem2D &S,
    real_t *Aflat,
    real_t *b)
{
    auto A = [&](int r, int c) -> real_t &
    {
        return Aflat[r * 5 + c];
    };

    const real_t rho_part =
        S.rho.lin.template get<MomentId2D::rho>() +
        S.rho.nonlin.template get<NonlinearId2D::uxux>() * r::zero + // placeholder
        S.rho.nonlin.template get<NonlinearId2D::uxuy>() * r::zero +
        S.rho.nonlin.template get<NonlinearId2D::uyuy>() * r::zero;

    const real_t uxI = S.incomings.template get<MomentId2D::ux>();
    const real_t uyI = S.incomings.template get<MomentId2D::uy>();
    const real_t mxxI = S.incomings.template get<MomentId2D::mxx>();
    const real_t mxyI = S.incomings.template get<MomentId2D::mxy>();
    const real_t myyI = S.incomings.template get<MomentId2D::myy>();

    // Linha ux
    A(0, 0) = S.ux.template get<MomentId2D::ux>() - S.rho.lin.template get<MomentId2D::ux>() * uxI;
    A(0, 1) = S.ux.template get<MomentId2D::uy>() - S.rho.lin.template get<MomentId2D::uy>() * uxI;
    A(0, 2) = S.ux.template get<MomentId2D::mxx>() - S.rho.lin.template get<MomentId2D::mxx>() * uxI;
    A(0, 3) = S.ux.template get<MomentId2D::mxy>() - S.rho.lin.template get<MomentId2D::mxy>() * uxI;
    A(0, 4) = S.ux.template get<MomentId2D::myy>() - S.rho.lin.template get<MomentId2D::myy>() * uxI;
    b[0] = uxI * rho_part - S.ux.template get<MomentId2D::rho>();

    // Linha uy
    A(1, 0) = S.uy.template get<MomentId2D::ux>() - S.rho.lin.template get<MomentId2D::ux>() * uyI;
    A(1, 1) = S.uy.template get<MomentId2D::uy>() - S.rho.lin.template get<MomentId2D::uy>() * uyI;
    A(1, 2) = S.uy.template get<MomentId2D::mxx>() - S.rho.lin.template get<MomentId2D::mxx>() * uyI;
    A(1, 3) = S.uy.template get<MomentId2D::mxy>() - S.rho.lin.template get<MomentId2D::mxy>() * uyI;
    A(1, 4) = S.uy.template get<MomentId2D::myy>() - S.rho.lin.template get<MomentId2D::myy>() * uyI;
    b[1] = uyI * rho_part - S.uy.template get<MomentId2D::rho>();

    // Linha mxx
    A(2, 0) = S.mxx.template get<MomentId2D::ux>() - S.rho.lin.template get<MomentId2D::ux>() * mxxI;
    A(2, 1) = S.mxx.template get<MomentId2D::uy>() - S.rho.lin.template get<MomentId2D::uy>() * mxxI;
    A(2, 2) = S.mxx.template get<MomentId2D::mxx>() - S.rho.lin.template get<MomentId2D::mxx>() * mxxI;
    A(2, 3) = S.mxx.template get<MomentId2D::mxy>() - S.rho.lin.template get<MomentId2D::mxy>() * mxxI;
    A(2, 4) = S.mxx.template get<MomentId2D::myy>() - S.rho.lin.template get<MomentId2D::myy>() * mxxI;
    b[2] = mxxI * rho_part - S.mxx.template get<MomentId2D::rho>();

    // Linha mxy
    A(3, 0) = S.mxy.template get<MomentId2D::ux>() - S.rho.lin.template get<MomentId2D::ux>() * mxyI;
    A(3, 1) = S.mxy.template get<MomentId2D::uy>() - S.rho.lin.template get<MomentId2D::uy>() * mxyI;
    A(3, 2) = S.mxy.template get<MomentId2D::mxx>() - S.rho.lin.template get<MomentId2D::mxx>() * mxyI;
    A(3, 3) = S.mxy.template get<MomentId2D::mxy>() - S.rho.lin.template get<MomentId2D::mxy>() * mxyI;
    A(3, 4) = S.mxy.template get<MomentId2D::myy>() - S.rho.lin.template get<MomentId2D::myy>() * mxyI;
    b[3] = mxyI * rho_part - S.mxy.template get<MomentId2D::rho>();

    // Linha myy
    A(4, 0) = S.myy.template get<MomentId2D::ux>() - S.rho.lin.template get<MomentId2D::ux>() * myyI;
    A(4, 1) = S.myy.template get<MomentId2D::uy>() - S.rho.lin.template get<MomentId2D::uy>() * myyI;
    A(4, 2) = S.myy.template get<MomentId2D::mxx>() - S.rho.lin.template get<MomentId2D::mxx>() * myyI;
    A(4, 3) = S.myy.template get<MomentId2D::mxy>() - S.rho.lin.template get<MomentId2D::mxy>() * myyI;
    A(4, 4) = S.myy.template get<MomentId2D::myy>() - S.rho.lin.template get<MomentId2D::myy>() * myyI;
    b[4] = myyI * rho_part - S.myy.template get<MomentId2D::rho>();
}