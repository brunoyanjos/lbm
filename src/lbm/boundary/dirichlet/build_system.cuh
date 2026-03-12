#pragma once

#include "../../../core/types.cuh"
#include "system.cuh"

__device__ __forceinline__ void build_dirichlet_system(
    const DirichletSystem &S,
    real_t ux,
    real_t uy,
    real_t *Aflat,
    real_t *b)
{
    auto A = [&](int r, int c) -> real_t &
    {
        return Aflat[r * 3 + c];
    };

    const real_t rho_part =
        S.rho.lin.template get<MomentId::rho>() +
        ux * S.rho.lin.template get<MomentId::ux>() +
        uy * S.rho.lin.template get<MomentId::uy>() +
        ux * ux * S.rho.nonlin.template get<NonlinearId::uxux>() +
        ux * uy * S.rho.nonlin.template get<NonlinearId::uxuy>() +
        uy * uy * S.rho.nonlin.template get<NonlinearId::uyuy>();

    const real_t mxxI = S.incomings.template get<MomentId::mxx>();
    const real_t mxyI = S.incomings.template get<MomentId::mxy>();
    const real_t myyI = S.incomings.template get<MomentId::myy>();

    A(0, 0) = S.mxx.template get<MomentId::mxx>() - S.rho.lin.template get<MomentId::mxx>() * mxxI;
    A(0, 1) = S.mxx.template get<MomentId::mxy>() - S.rho.lin.template get<MomentId::mxy>() * mxxI;
    A(0, 2) = S.mxx.template get<MomentId::myy>() - S.rho.lin.template get<MomentId::myy>() * mxxI;
    b[0] = mxxI * rho_part - (S.mxx.template get<MomentId::rho>() +
                              ux * S.mxx.template get<MomentId::ux>() + uy * S.mxx.template get<MomentId::uy>());

    A(1, 0) = S.mxy.template get<MomentId::mxx>() - S.rho.lin.template get<MomentId::mxx>() * mxyI;
    A(1, 1) = S.mxy.template get<MomentId::mxy>() - S.rho.lin.template get<MomentId::mxy>() * mxyI;
    A(1, 2) = S.mxy.template get<MomentId::myy>() - S.rho.lin.template get<MomentId::myy>() * mxyI;
    b[1] = mxyI * rho_part - (S.mxy.template get<MomentId::rho>() +
                              ux * S.mxy.template get<MomentId::ux>() + uy * S.mxy.template get<MomentId::uy>());

    A(2, 0) = S.myy.template get<MomentId::mxx>() - S.rho.lin.template get<MomentId::mxx>() * myyI;
    A(2, 1) = S.myy.template get<MomentId::mxy>() - S.rho.lin.template get<MomentId::mxy>() * myyI;
    A(2, 2) = S.myy.template get<MomentId::myy>() - S.rho.lin.template get<MomentId::myy>() * myyI;
    b[2] = myyI * rho_part - (S.myy.template get<MomentId::rho>() +
                              ux * S.myy.template get<MomentId::ux>() + uy * S.myy.template get<MomentId::uy>());
}