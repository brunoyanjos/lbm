#pragma once

#include "../../../core/types.cuh"
#include "system.cuh"
#include "newton_coefficients.cuh"

__device__ __forceinline__ void build_fluid_newton_coefficients(
    const FluidSystem &S,
    FluidNewtonCoefficients &C)
{
    auto A = [&](int r, int c) -> real_t &
    { return C.coeff(r, c); };

    const real_t rho_rho = S.rho.lin.template get<MomentId::rho>();
    const real_t rho_ux = S.rho.lin.template get<MomentId::ux>();
    const real_t rho_uy = S.rho.lin.template get<MomentId::uy>();
    const real_t rho_mxx = S.rho.lin.template get<MomentId::mxx>();
    const real_t rho_mxy = S.rho.lin.template get<MomentId::mxy>();
    const real_t rho_myy = S.rho.lin.template get<MomentId::myy>();

    const real_t rho_uxux = S.rho.nonlin.template get<NonlinearId::uxux>();
    const real_t rho_uxuy = S.rho.nonlin.template get<NonlinearId::uxuy>();
    const real_t rho_uyuy = S.rho.nonlin.template get<NonlinearId::uyuy>();

    const real_t ux_I = S.incomings.template get<MomentId::ux>();
    const real_t uy_I = S.incomings.template get<MomentId::uy>();
    const real_t mxx_I = S.incomings.template get<MomentId::mxx>();
    const real_t mxy_I = S.incomings.template get<MomentId::mxy>();
    const real_t myy_I = S.incomings.template get<MomentId::myy>();

    // ux_I equation
    A(0, 0) = S.ux.template get<MomentId::ux>() - rho_ux * ux_I;
    A(0, 1) = S.ux.template get<MomentId::uy>() - rho_uy * ux_I;
    A(0, 2) = -rho_uxux * ux_I;
    A(0, 3) = -rho_uxuy * ux_I;
    A(0, 4) = -rho_uyuy * ux_I;
    A(0, 5) = S.ux.template get<MomentId::mxx>() - rho_mxx * ux_I;
    A(0, 6) = S.ux.template get<MomentId::mxy>() - rho_mxy * ux_I;
    A(0, 7) = S.ux.template get<MomentId::myy>() - rho_myy * ux_I;
    C.b[0] = ux_I * rho_rho - S.ux.template get<MomentId::rho>();

    // uy_I equation
    A(1, 0) = S.uy.template get<MomentId::ux>() - rho_ux * uy_I;
    A(1, 1) = S.uy.template get<MomentId::uy>() - rho_uy * uy_I;
    A(1, 2) = -rho_uxux * uy_I;
    A(1, 3) = -rho_uxuy * uy_I;
    A(1, 4) = -rho_uyuy * uy_I;
    A(1, 5) = S.uy.template get<MomentId::mxx>() - rho_mxx * uy_I;
    A(1, 6) = S.uy.template get<MomentId::mxy>() - rho_mxy * uy_I;
    A(1, 7) = S.uy.template get<MomentId::myy>() - rho_myy * uy_I;
    C.b[1] = uy_I * rho_rho - S.uy.template get<MomentId::rho>();

    // mxx_I equation
    A(2, 0) = S.mxx.template get<MomentId::ux>() - rho_ux * mxx_I;
    A(2, 1) = S.mxx.template get<MomentId::uy>() - rho_uy * mxx_I;
    A(2, 2) = -rho_uxux * mxx_I;
    A(2, 3) = -rho_uxuy * mxx_I;
    A(2, 4) = -rho_uyuy * mxx_I;
    A(2, 5) = S.mxx.template get<MomentId::mxx>() - rho_mxx * mxx_I;
    A(2, 6) = S.mxx.template get<MomentId::mxy>() - rho_mxy * mxx_I;
    A(2, 7) = S.mxx.template get<MomentId::myy>() - rho_myy * mxx_I;
    C.b[2] = mxx_I * rho_rho - S.mxx.template get<MomentId::rho>();

    // mxy_I equation
    A(3, 0) = S.mxy.template get<MomentId::ux>() - rho_ux * mxy_I;
    A(3, 1) = S.mxy.template get<MomentId::uy>() - rho_uy * mxy_I;
    A(3, 2) = -rho_uxux * mxy_I;
    A(3, 3) = -rho_uxuy * mxy_I;
    A(3, 4) = -rho_uyuy * mxy_I;
    A(3, 5) = S.mxy.template get<MomentId::mxx>() - rho_mxx * mxy_I;
    A(3, 6) = S.mxy.template get<MomentId::mxy>() - rho_mxy * mxy_I;
    A(3, 7) = S.mxy.template get<MomentId::myy>() - rho_myy * mxy_I;
    C.b[3] = mxy_I * rho_rho - S.mxy.template get<MomentId::rho>();

    // myy_I equation
    A(4, 0) = S.myy.template get<MomentId::ux>() - rho_ux * myy_I;
    A(4, 1) = S.myy.template get<MomentId::uy>() - rho_uy * myy_I;
    A(4, 2) = -rho_uxux * myy_I;
    A(4, 3) = -rho_uxuy * myy_I;
    A(4, 4) = -rho_uyuy * myy_I;
    A(4, 5) = S.myy.template get<MomentId::mxx>() - rho_mxx * myy_I;
    A(4, 6) = S.myy.template get<MomentId::mxy>() - rho_mxy * myy_I;
    A(4, 7) = S.myy.template get<MomentId::myy>() - rho_myy * myy_I;
    C.b[4] = myy_I * rho_rho - S.myy.template get<MomentId::rho>();
}