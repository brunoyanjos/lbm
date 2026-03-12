#pragma once

#include "../../../core/types.cuh"
#include "../common/moment_ids.cuh"
#include "../common/nonlinear_ids.cuh"
#include "system.cuh"

__device__ __forceinline__ void build_fluid_system(
    const FluidSystem &S,
    real_t *Aflat,
    real_t *b)
{
    auto A = [&](int r, int c) -> real_t &
    {
        return Aflat[r * 5 + c];
    };

    const real_t rho_part =
        S.rho.lin.template get<MomentId::rho>() +
        S.rho.nonlin.template get<NonlinearId::uxux>() * r::zero + // placeholder
        S.rho.nonlin.template get<NonlinearId::uxuy>() * r::zero +
        S.rho.nonlin.template get<NonlinearId::uyuy>() * r::zero;

    const real_t uxI = S.incomings.template get<MomentId::ux>();
    const real_t uyI = S.incomings.template get<MomentId::uy>();
    const real_t mxxI = S.incomings.template get<MomentId::mxx>();
    const real_t mxyI = S.incomings.template get<MomentId::mxy>();
    const real_t myyI = S.incomings.template get<MomentId::myy>();

    // Linha ux
    A(0, 0) = S.ux.template get<MomentId::ux>() - S.rho.lin.template get<MomentId::ux>() * uxI;
    A(0, 1) = S.ux.template get<MomentId::uy>() - S.rho.lin.template get<MomentId::uy>() * uxI;
    A(0, 2) = S.ux.template get<MomentId::mxx>() - S.rho.lin.template get<MomentId::mxx>() * uxI;
    A(0, 3) = S.ux.template get<MomentId::mxy>() - S.rho.lin.template get<MomentId::mxy>() * uxI;
    A(0, 4) = S.ux.template get<MomentId::myy>() - S.rho.lin.template get<MomentId::myy>() * uxI;
    b[0] = uxI * rho_part - S.ux.template get<MomentId::rho>();

    // Linha uy
    A(1, 0) = S.uy.template get<MomentId::ux>() - S.rho.lin.template get<MomentId::ux>() * uyI;
    A(1, 1) = S.uy.template get<MomentId::uy>() - S.rho.lin.template get<MomentId::uy>() * uyI;
    A(1, 2) = S.uy.template get<MomentId::mxx>() - S.rho.lin.template get<MomentId::mxx>() * uyI;
    A(1, 3) = S.uy.template get<MomentId::mxy>() - S.rho.lin.template get<MomentId::mxy>() * uyI;
    A(1, 4) = S.uy.template get<MomentId::myy>() - S.rho.lin.template get<MomentId::myy>() * uyI;
    b[1] = uyI * rho_part - S.uy.template get<MomentId::rho>();

    // Linha mxx
    A(2, 0) = S.mxx.template get<MomentId::ux>() - S.rho.lin.template get<MomentId::ux>() * mxxI;
    A(2, 1) = S.mxx.template get<MomentId::uy>() - S.rho.lin.template get<MomentId::uy>() * mxxI;
    A(2, 2) = S.mxx.template get<MomentId::mxx>() - S.rho.lin.template get<MomentId::mxx>() * mxxI;
    A(2, 3) = S.mxx.template get<MomentId::mxy>() - S.rho.lin.template get<MomentId::mxy>() * mxxI;
    A(2, 4) = S.mxx.template get<MomentId::myy>() - S.rho.lin.template get<MomentId::myy>() * mxxI;
    b[2] = mxxI * rho_part - S.mxx.template get<MomentId::rho>();

    // Linha mxy
    A(3, 0) = S.mxy.template get<MomentId::ux>() - S.rho.lin.template get<MomentId::ux>() * mxyI;
    A(3, 1) = S.mxy.template get<MomentId::uy>() - S.rho.lin.template get<MomentId::uy>() * mxyI;
    A(3, 2) = S.mxy.template get<MomentId::mxx>() - S.rho.lin.template get<MomentId::mxx>() * mxyI;
    A(3, 3) = S.mxy.template get<MomentId::mxy>() - S.rho.lin.template get<MomentId::mxy>() * mxyI;
    A(3, 4) = S.mxy.template get<MomentId::myy>() - S.rho.lin.template get<MomentId::myy>() * mxyI;
    b[3] = mxyI * rho_part - S.mxy.template get<MomentId::rho>();

    // Linha myy
    A(4, 0) = S.myy.template get<MomentId::ux>() - S.rho.lin.template get<MomentId::ux>() * myyI;
    A(4, 1) = S.myy.template get<MomentId::uy>() - S.rho.lin.template get<MomentId::uy>() * myyI;
    A(4, 2) = S.myy.template get<MomentId::mxx>() - S.rho.lin.template get<MomentId::mxx>() * myyI;
    A(4, 3) = S.myy.template get<MomentId::mxy>() - S.rho.lin.template get<MomentId::mxy>() * myyI;
    A(4, 4) = S.myy.template get<MomentId::myy>() - S.rho.lin.template get<MomentId::myy>() * myyI;
    b[4] = myyI * rho_part - S.myy.template get<MomentId::rho>();
}