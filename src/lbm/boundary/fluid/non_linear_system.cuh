#pragma once

#include "../../../core/types.cuh"
#include "../../../core/math_utils.cuh"
#include "../../../core/linear_solver.cuh"
#include "system_variables.cuh"

struct FluidNewtonCoefficients2D
{
    real_t A[5 * 8]{};
    real_t b[5]{};

    __device__ __forceinline__ real_t &coeff(int r, int c)
    {
        return A[r * 8 + c];
    }

    __device__ __forceinline__ const real_t &coeff(int r, int c) const
    {
        return A[r * 8 + c];
    }
};

__device__ __forceinline__ void build_fluid_newton_coefficients(
    const FluidSystem2D &S,
    FluidNewtonCoefficients2D &C)
{
    auto A = [&](int r, int c) -> real_t &
    { return C.coeff(r, c); };

    // ux_I equation
    A(0, 0) = S.ux_ux - S.rho_ux * S.ux_I;   // ux
    A(0, 1) = S.ux_uy - S.rho_uy * S.ux_I;   // uy
    A(0, 2) = -S.rho_uxux * S.ux_I;          // ux ux
    A(0, 3) = -S.rho_uxuy * S.ux_I;          // ux uy
    A(0, 4) = -S.rho_uyuy * S.ux_I;          // uy uy
    A(0, 5) = S.ux_mxx - S.rho_mxx * S.ux_I; // mxx
    A(0, 6) = S.ux_mxy - S.rho_mxy * S.ux_I; // mxy
    A(0, 7) = S.ux_myy - S.rho_myy * S.ux_I; // myy
    C.b[0] = S.ux_I * S.rho_rho - S.ux_rho;

    // uy_I equation
    A(1, 0) = S.uy_ux - S.rho_ux * S.uy_I;   // ux
    A(1, 1) = S.uy_uy - S.rho_uy * S.uy_I;   // uy
    A(1, 2) = -S.rho_uxux * S.uy_I;          // ux ux
    A(1, 3) = -S.rho_uxuy * S.uy_I;          // ux uy
    A(1, 4) = -S.rho_uyuy * S.uy_I;          // uy uy
    A(1, 5) = S.uy_mxx - S.rho_mxx * S.uy_I; // mxx
    A(1, 6) = S.uy_mxy - S.rho_mxy * S.uy_I; // mxy
    A(1, 7) = S.uy_myy - S.rho_myy * S.uy_I; // myy
    C.b[1] = S.uy_I * S.rho_rho - S.uy_rho;

    // mxx_I equation
    A(2, 0) = S.mxx_ux - S.rho_ux * S.mxx_I;   // ux
    A(2, 1) = S.mxx_uy - S.rho_uy * S.mxx_I;   // uy
    A(2, 2) = -S.rho_uxux * S.mxx_I;           // ux ux
    A(2, 3) = -S.rho_uxuy * S.mxx_I;           // ux uy
    A(2, 4) = -S.rho_uyuy * S.mxx_I;           // uy uy
    A(2, 5) = S.mxx_mxx - S.rho_mxx * S.mxx_I; // mxx
    A(2, 6) = S.mxx_mxy - S.rho_mxy * S.mxx_I; // mxy
    A(2, 7) = S.mxx_myy - S.rho_myy * S.mxx_I; // myy
    C.b[2] = S.mxx_I * S.rho_rho - S.mxx_rho;

    // mxy_I equation
    A(3, 0) = S.mxy_ux - S.rho_ux * S.mxy_I;   // ux
    A(3, 1) = S.mxy_uy - S.rho_uy * S.mxy_I;   // uy
    A(3, 2) = -S.rho_uxux * S.mxy_I;           // ux ux
    A(3, 3) = -S.rho_uxuy * S.mxy_I;           // ux uy
    A(3, 4) = -S.rho_uyuy * S.mxy_I;           // uy uy
    A(3, 5) = S.mxy_mxx - S.rho_mxx * S.mxy_I; // mxx
    A(3, 6) = S.mxy_mxy - S.rho_mxy * S.mxy_I; // mxy
    A(3, 7) = S.mxy_myy - S.rho_myy * S.mxy_I; // myy
    C.b[3] = S.mxy_I * S.rho_rho - S.mxy_rho;

    // myy_I equation
    A(4, 0) = S.myy_ux - S.rho_ux * S.myy_I;   // ux
    A(4, 1) = S.myy_uy - S.rho_uy * S.myy_I;   // uy
    A(4, 2) = -S.rho_uxux * S.myy_I;           // ux ux
    A(4, 3) = -S.rho_uxuy * S.myy_I;           // ux uy
    A(4, 4) = -S.rho_uyuy * S.myy_I;           // uy uy
    A(4, 5) = S.myy_mxx - S.rho_mxx * S.myy_I; // mxx
    A(4, 6) = S.myy_mxy - S.rho_mxy * S.myy_I; // mxy
    A(4, 7) = S.myy_myy - S.rho_myy * S.myy_I; // myy
    C.b[4] = S.myy_I * S.rho_rho - S.myy_rho;
}

__device__ __forceinline__ real_t rel_step(real_t x_new, real_t x_old)
{
    const real_t eps = real_t(1e-12);
    const real_t denom = fmax(r_abs(x_new), eps);
    return r_abs(x_new - x_old) / denom;
}

__device__ __forceinline__ real_t eval_fluid_row(
    const real_t *__restrict__ row,
    real_t b_row,
    real_t ux,
    real_t uy,
    real_t mxx,
    real_t mxy,
    real_t myy)
{
    return row[0] * ux + row[1] * uy +
           row[2] * ux * ux + row[3] * ux * uy + row[4] * uy * uy +
           row[5] * mxx + row[6] * mxy + row[7] * myy -
           b_row;
}

__device__ __forceinline__ void build_fluid_jacobian(
    const FluidNewtonCoefficients2D &C,
    real_t ux,
    real_t uy,
    real_t *__restrict__ J)
{
    auto U = [&](int r, int c) -> real_t &
    { return J[r * 5 + c]; };

#pragma unroll
    for (int i = 0; i < 5; ++i)
    {
        U(i, 0) = C.coeff(i, 0) + r::two * C.coeff(i, 2) * ux + C.coeff(i, 3) * uy; // df/dux
        U(i, 1) = C.coeff(i, 1) + C.coeff(i, 3) * ux + r::two * C.coeff(i, 4) * uy; // df/duy
        U(i, 2) = C.coeff(i, 5);                                                    // df/dmxx
        U(i, 3) = C.coeff(i, 6);                                                    // df/dmxy
        U(i, 4) = C.coeff(i, 7);                                                    // df/dmyy
    }
}

__device__ __forceinline__ void build_fluid_rhs(
    const FluidNewtonCoefficients2D &C,
    real_t ux,
    real_t uy,
    real_t mxx,
    real_t mxy,
    real_t myy,
    real_t *__restrict__ rhs)
{
#pragma unroll
    for (int i = 0; i < 5; ++i)
        rhs[i] = -eval_fluid_row(&C.A[i * 8], C.b[i], ux, uy, mxx, mxy, myy);
}

__device__ __forceinline__ void solve_fluid_newton(
    const FluidNewtonCoefficients2D &C,
    real_t &ux,
    real_t &uy,
    real_t &mxx,
    real_t &mxy,
    real_t &myy)
{
    real_t J[5 * 5]{};
    real_t rhs[5]{};
    real_t delta[5]{};

    real_t error = r::one;
    int it = 0;
    constexpr int it_max = 50;
    constexpr real_t tol = real_t(1e-6);

    while (error > tol && it++ < it_max)
    {
        const real_t ux_old = ux;
        const real_t uy_old = uy;
        const real_t mxx_old = mxx;
        const real_t mxy_old = mxy;
        const real_t myy_old = myy;

        build_fluid_jacobian(C, ux, uy, J);
        build_fluid_rhs(C, ux, uy, mxx, mxy, myy, rhs);

        gaussianElimination<5>(J, rhs, delta);

        ux = ux_old + delta[0];
        uy = uy_old + delta[1];
        mxx = mxx_old + delta[2];
        mxy = mxy_old + delta[3];
        myy = myy_old + delta[4];

        error = rel_step(ux, ux_old);
        error = fmax(error, rel_step(uy, uy_old));
        error = fmax(error, rel_step(mxx, mxx_old));
        error = fmax(error, rel_step(mxy, mxy_old));
        error = fmax(error, rel_step(myy, myy_old));
    }
}