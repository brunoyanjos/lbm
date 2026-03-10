#pragma once

#include "../../core/types.cuh"
#include "../../core/linear_solver.cuh"
#include "../domain/mask_utils.cuh"
#include "../../geometries/active_geometry.cuh"
#include "../../core/math_utils.cuh"
#include <cstdint>

#define STRONG_CONSERVATION_UX 0
#define STRONG_CONSERVATION_UY 0

#define EQUATION_ON_MXY 1

#define MXX_INCOMPRESSIBLE 1
#define MYY_INCOMPRESSIBLE 1

#define MXX_MINUS_MYY_INCOMPRESSIBLE 0

__device__ inline void evaluate_fluid_boundary(
    real_t *__restrict__ pop,
    uint32_t valid_mask,
    real_t &rho,
    real_t &ux, real_t &uy,
    real_t &mxx, real_t &mxy, real_t &myy,
    int x, int y)
{
    const uint32_t outgoing_mask = valid_mask;
    const uint32_t incoming_mask = mask_opp(valid_mask);

    real_t A_coeff[40]{};
    real_t b_coeff[5]{};

    real_t rho_I = r::zero;

    real_t ux_I = r::zero;
    real_t uy_I = r::zero;

    real_t mxx_I = r::zero;
    real_t mxy_I = r::zero;
    real_t myy_I = r::zero;

    real_t A_rho[9]{};

    // real_t A_ux[6]{};
    // real_t A_uy[6]{};

    real_t A_ux[9]{};
    real_t A_uy[9]{};

    real_t A_mxx[6]{};
    real_t A_mxy[6]{};
    real_t A_myy[6]{};

    real_t c, s, r;
    polar_unit_vectors(x, y, x, y, c, s, r);

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        real_t cx, cy, Hxx, Hxy, Hyy;
        Stencil::basis2polar(i, c, s, cx, cy, Hxx, Hxy, Hyy);

        const real_t factors[6] = {
            Stencil::w(i),
            Stencil::as2 * Stencil::w(i) * cx,
            Stencil::as2 * Stencil::w(i) * cy,
            Stencil::as4 * r::half * Stencil::w(i) * Hxx,
            Stencil::as4 * Stencil::w(i) * Hxy,
            Stencil::as4 * r::half * Stencil::w(i) * Hyy};

        if (dir_valid(incoming_mask, i))
        {
            rho_I += pop[i];
            ux_I += pop[i] * cx;
            uy_I += pop[i] * cy;
            mxx_I += pop[i] * Hxx;
            mxy_I += pop[i] * Hxy;
            myy_I += pop[i] * Hyy;

#pragma unroll
            for (int j = 0; j < 6; ++j)
            {
#if !STRONG_CONSERVATION_UX
                A_ux[j] += factors[j] * cx;
#endif

#if !STRONG_CONSERVATION_UY
                A_uy[j] += factors[j] * cy;
#endif

                A_mxx[j] += factors[j] * Hxx;
                A_mxy[j] += factors[j] * Hxy;
                A_myy[j] += factors[j] * Hyy;
            }
        }

        if (dir_valid(outgoing_mask, i))
        {
            A_rho[0] += factors[0];
            A_rho[1] += factors[1];
            A_rho[2] += factors[2];
            A_rho[3] += Geometry::OMEGA * factors[3];
            A_rho[4] += Geometry::OMEGA * factors[4];
            A_rho[5] += Geometry::OMEGA * factors[5];
            A_rho[6] += (r::one - Geometry::OMEGA) * factors[3];
            A_rho[7] += (r::one - Geometry::OMEGA) * factors[4];
            A_rho[8] += (r::one - Geometry::OMEGA) * factors[5];

#if STRONG_CONSERVATION_UX
            A_ux[0] += factors[0] * cx;
            A_ux[1] += factors[1] * cx;
            A_ux[2] += factors[2] * cx;
            A_ux[3] += OMEGA * factors[3] * cx;
            A_ux[4] += OMEGA * factors[4] * cx;
            A_ux[5] += OMEGA * factors[5] * cx;
            A_ux[6] += (r::one - OMEGA) * factors[3] * cx;
            A_ux[7] += (r::one - OMEGA) * factors[4] * cx;
            A_ux[8] += (r::one - OMEGA) * factors[5] * cx;
#endif

#if STRONG_CONSERVATION_UY
            A_uy[0] += factors[0] * cy;
            A_uy[1] += factors[1] * cy;
            A_uy[2] += factors[2] * cy;
            A_uy[3] += OMEGA * factors[3] * cy;
            A_uy[4] += OMEGA * factors[4] * cy;
            A_uy[5] += OMEGA * factors[5] * cy;
            A_uy[6] += (r::one - OMEGA) * factors[3] * cy;
            A_uy[7] += (r::one - OMEGA) * factors[4] * cy;
            A_uy[8] += (r::one - OMEGA) * factors[5] * cy;
#endif
        }
    }

    const real_t inv_rho_I = real_t(1) / rho_I;

    ux_I *= inv_rho_I;
    uy_I *= inv_rho_I;
    mxx_I *= inv_rho_I;
    mxy_I *= inv_rho_I;
    myy_I *= inv_rho_I;

    auto A = [&](int r, int c) -> real_t &
    { return A_coeff[r * 8 + c]; };

    // uxI equation
#if EQUATION_ON_MXY
#if !STRONG_CONSERVATION_UX
    A(0, 0) = A_ux[1] - A_rho[1] * ux_I; // ux
    A(0, 1) = A_ux[2] - A_rho[2] * ux_I; // uy
    A(0, 2) = -A_rho[3] * ux_I;          // ux ux
    A(0, 3) = -A_rho[4] * ux_I;          // ux uy
    A(0, 4) = -A_rho[5] * ux_I;          // uy uy
    A(0, 5) = A_ux[3] - A_rho[6] * ux_I; // mxx
    A(0, 6) = A_ux[4] - A_rho[7] * ux_I; // mxy
    A(0, 7) = A_ux[5] - A_rho[8] * ux_I; // myy

    b_coeff[0] = ux_I * A_rho[0] - A_ux[0];
#else
    A(0, 0) = A_ux[1] + A_rho[1] * ux_I; // ux
    A(0, 1) = A_ux[2] + A_rho[2] * ux_I; // uy
    A(0, 2) = A_ux[3] + A_rho[3] * ux_I; // ux ux
    A(0, 3) = A_ux[4] + A_rho[4] * ux_I; // ux uy
    A(0, 4) = A_ux[5] + A_rho[5] * ux_I; // uy uy
    A(0, 5) = A_ux[6] + A_rho[6] * ux_I; // mxx
    A(0, 6) = A_ux[7] + A_rho[7] * ux_I; // mxy
    A(0, 7) = A_ux[8] + A_rho[8] * ux_I; // myy

    b_coeff[0] = -ux_I * A_rho[0] - A_ux[0];
#endif
#else
    // const real_t delta = r_cast(0.1);

    // if (y < NY / 2)
    // {
    //     A(0, 0) = 1;                             // ux
    //     A(0, 1) = r::zero;                       // uy
    //     A(0, 2) = r::zero;                       // ux ux
    //     A(0, 3) = -Stencil::as2 * OMEGA * delta; // ux uy
    //     A(0, 4) = r::zero;                       // uy uy
    //     A(0, 5) = r::zero;                       // mxx
    //     A(0, 6) = Stencil::as2 * OMEGA * delta;  // mxy
    //     A(0, 7) = r::zero;                       // myy

    //     b_coeff[1] = r::zero;
    // }
    // else
    // {
    //     A(0, 0) = 1;                             // ux
    //     A(0, 1) = r::zero;                       // uy
    //     A(0, 2) = r::zero;                       // ux ux
    //     A(0, 3) = Stencil::as2 * OMEGA * delta;  // ux uy
    //     A(0, 4) = r::zero;                       // uy uy
    //     A(0, 5) = r::zero;                       // mxx
    //     A(0, 6) = -Stencil::as2 * OMEGA * delta; // mxy
    //     A(0, 7) = r::zero;                       // myy

    //     b_coeff[1] = U_MAX;
    // }
#endif

    // uyI equation
#if !STRONG_CONSERVATION_UY
    A(1, 0) = A_uy[1] - A_rho[1] * uy_I; // ux
    A(1, 1) = A_uy[2] - A_rho[2] * uy_I; // uy
    A(1, 2) = -A_rho[3] * uy_I;          // ux ux
    A(1, 3) = -A_rho[4] * uy_I;          // ux uy
    A(1, 4) = -A_rho[5] * uy_I;          // uy uy
    A(1, 5) = A_uy[3] - A_rho[6] * uy_I; // mxx
    A(1, 6) = A_uy[4] - A_rho[7] * uy_I; // mxy
    A(1, 7) = A_uy[5] - A_rho[8] * uy_I; // myy

    b_coeff[1] = uy_I * A_rho[0] - A_uy[0];
#else
    A(1, 0) = A_uy[1] + A_rho[1] * uy_I; // ux
    A(1, 1) = A_uy[2] + A_rho[2] * uy_I; // uy
    A(1, 2) = A_uy[3] + A_rho[3] * uy_I; // ux ux
    A(1, 3) = A_uy[4] + A_rho[4] * uy_I; // ux uy
    A(1, 4) = A_uy[5] + A_rho[5] * uy_I; // uy uy
    A(1, 5) = A_uy[6] + A_rho[6] * uy_I; // mxx
    A(1, 6) = A_uy[7] + A_rho[7] * uy_I; // mxy
    A(1, 7) = A_uy[8] + A_rho[8] * uy_I; // myy

    b_coeff[1] = -uy_I * A_rho[0] - A_uy[0];
#endif

// mxxI equation
#if !MXX_INCOMPRESSIBLE
    A(2, 0) = A_mxx[1] - A_rho[1] * mxx_I; // ux
    A(2, 1) = A_mxx[2] - A_rho[2] * mxx_I; // uy
    A(2, 2) = -A_rho[3] * mxx_I;           // ux ux
    A(2, 3) = -A_rho[4] * mxx_I;           // ux uy
    A(2, 4) = -A_rho[5] * mxx_I;           // uy uy
    A(2, 5) = A_mxx[3] - A_rho[6] * mxx_I; // mxx
    A(2, 6) = A_mxx[4] - A_rho[7] * mxx_I; // mxy
    A(2, 7) = A_mxx[5] - A_rho[8] * mxx_I; // myy

    b_coeff[2] = mxx_I * A_rho[0] - A_mxx[0];
#else
#if MXX_MINUS_MYY_INCOMPRESSIBLE
    A(2, 0) = A_mxx[1] - A_rho[1] * mxx_I - (A_myy[1] - A_rho[1] * myy_I); // ux
    A(2, 1) = A_mxx[2] - A_rho[2] * mxx_I - (A_myy[2] - A_rho[2] * myy_I); // uy
    A(2, 2) = A_rho[3] * myy_I - A_rho[3] * mxx_I;                         // ux ux
    A(2, 3) = A_rho[4] * myy_I - A_rho[4] * mxx_I;                         // ux uy
    A(2, 4) = A_rho[5] * myy_I - A_rho[5] * mxx_I;                         // uy uy
    A(2, 5) = A_mxx[3] - A_rho[6] * mxx_I - (A_myy[3] - A_rho[6] * myy_I); // mxx
    A(2, 6) = A_mxx[4] - A_rho[7] * mxx_I - (A_myy[4] - A_rho[7] * myy_I); // mxy
    A(2, 7) = A_mxx[5] - A_rho[8] * mxx_I - (A_myy[5] - A_rho[8] * myy_I); // myy

    b_coeff[2] = mxx_I * A_rho[0] - A_mxx[0] - (myy_I * A_rho[0] - A_myy[0]);
#else
    A(2, 0) = r::zero; // ux
    A(2, 1) = r::zero; // uy
    A(2, 2) = -1;      // ux ux
    A(2, 3) = r::zero; // ux uy
    A(2, 4) = r::zero; // uy uy
    A(2, 5) = 1;       // mxx
    A(2, 6) = r::zero; // mxy
    A(2, 7) = r::zero; // myy

    b_coeff[2] = r::zero;
#endif
#endif

    // mxyI equation
#if EQUATION_ON_MXY
    // if (r < (Geometry::R_IN + Geometry::R_OUT) * r::half)
    // {
    //     A(3, 0) = r::zero;                            // ux
    //     A(3, 1) = 1;                                  // uy
    //     A(3, 2) = r::zero;                            // ux ux
    //     A(3, 3) = -Stencil::as2 * Geometry::OMEGA * (r - R_IN); // ux uy
    //     A(3, 4) = r::zero;                            // uy uy
    //     A(3, 5) = r::zero;                            // mxx
    //     A(3, 6) = Stencil::as2 * Geometry::OMEGA * (r - R_IN);  // mxy
    //     A(3, 7) = r::zero;                            // myy

    //     b_coeff[3] = Geometry::U_MAX;
    // }
    // else
    // {
    //     A(3, 0) = r::zero;                             // ux
    //     A(3, 1) = 1;                                   // uy
    //     A(3, 2) = r::zero;                             // ux ux
    //     A(3, 3) = Stencil::as2 * Geometry::OMEGA * (Geometry::R_OUT - r);  // ux uy
    //     A(3, 4) = r::zero;                             // uy uy
    //     A(3, 5) = r::zero;                             // mxx
    //     A(3, 6) = -Stencil::as2 * Geometry::OMEGA * (Geometry::R_OUT - r); // mxy
    //     A(3, 7) = r::zero;                             // myy

    //     b_coeff[3] = r::zero;
    // }
#else
    A(3, 0) = A_mxy[1] - A_rho[1] * mxy_I; // ux
    A(3, 1) = A_mxy[2] - A_rho[2] * mxy_I; // uy
    A(3, 2) = -A_rho[3] * mxy_I;           // ux ux
    A(3, 3) = -A_rho[4] * mxy_I;           // ux uy
    A(3, 4) = -A_rho[5] * mxy_I;           // uy uy
    A(3, 5) = A_mxy[3] - A_rho[6] * mxy_I; // mxx
    A(3, 6) = A_mxy[4] - A_rho[7] * mxy_I; // mxy
    A(3, 7) = A_mxy[5] - A_rho[8] * mxy_I; // myy

    b_coeff[3] = mxy_I * A_rho[0] - A_mxy[0];
#endif

// myyI equation
#if !MYY_INCOMPRESSIBLE
    A(4, 0) = A_myy[1] - A_rho[1] * myy_I; // ux
    A(4, 1) = A_myy[2] - A_rho[2] * myy_I; // uy
    A(4, 2) = -A_rho[3] * myy_I;           // ux ux
    A(4, 3) = -A_rho[4] * myy_I;           // ux uy
    A(4, 4) = -A_rho[5] * myy_I;           // uy uy
    A(4, 5) = A_myy[3] - A_rho[6] * myy_I; // mxx
    A(4, 6) = A_myy[4] - A_rho[7] * myy_I; // mxy
    A(4, 7) = A_myy[5] - A_rho[8] * myy_I; // myy

    b_coeff[4] = myy_I * A_rho[0] - A_myy[0];
#else
#if MXX_MINUS_MYY_INCOMPRESSIBLE
    A(4, 0) = r::zero; // ux
    A(4, 1) = r::zero; // uy
    A(4, 2) = -1;      // ux ux
    A(4, 3) = r::zero; // ux uy
    A(4, 4) = -1;      // uy uy
    A(4, 5) = 1;       // mxx
    A(4, 6) = r::zero; // mxy
    A(4, 7) = 1;       // myy

    b_coeff[4] = r::zero;
#else
    A(4, 0) = r::zero; // ux
    A(4, 1) = r::zero; // uy
    A(4, 2) = r::zero; // ux ux
    A(4, 3) = r::zero; // ux uy
    A(4, 4) = -1;      // uy uy
    A(4, 5) = r::zero; // mxx
    A(4, 6) = r::zero; // mxy
    A(4, 7) = 1;       // myy

    b_coeff[4] = r::zero;
#endif
#endif

    real_t A_gauss[5 * 5]{};
    real_t b_gauss[5]{};
    real_t moms[5]{};

    auto U = [&](int r, int c) -> real_t &
    { return A_gauss[r * 5 + c]; };

    real_t error = real_t(1);
    int it = 0;
    const int it_max = 50;

    while (error > real_t(1e-8) && it++ < it_max)
    {
        const real_t ux_old = ux;
        const real_t uy_old = uy;
        const real_t mxx_old = mxx;
        const real_t mxy_old = mxy;
        const real_t myy_old = myy;

        for (int i = 0; i < 5; ++i)
        {
            U(i, 0) = A(i, 0) + real_t(2) * A(i, 2) * ux + A(i, 3) * uy; // df_dux
            U(i, 1) = A(i, 1) + A(i, 3) * ux + real_t(2) * A(i, 4) * uy; // df_duy
            U(i, 2) = A(i, 5);                                           // df_dmxx
            U(i, 3) = A(i, 6);                                           // df_dmxy
            U(i, 4) = A(i, 7);                                           // df_dmyy

            b_gauss[i] = -eval_row(&A_coeff[i * 8], b_coeff[i], ux, uy, mxx, mxy, myy);
        }

        gaussianElimination<5>(A_gauss, b_gauss, moms);

        ux = ux_old + moms[0];
        uy = uy_old + moms[1];
        mxx = mxx_old + moms[2];
        mxy = mxy_old + moms[3];
        myy = myy_old + moms[4];

        // Erro = maior passo relativo
        error = rel_step(ux, ux_old);
        error = fmax(error, rel_step(uy, uy_old));
        error = fmax(error, rel_step(mxx, mxx_old));
        error = fmax(error, rel_step(mxy, mxy_old));
        error = fmax(error, rel_step(myy, myy_old));
    }

    real_t rho_denom = A_rho[0] +
                       A_rho[1] * ux +
                       A_rho[2] * uy +
                       A_rho[3] * ux * ux +
                       A_rho[4] * ux * uy +
                       A_rho[5] * uy * uy +
                       A_rho[6] * mxx +
                       A_rho[7] * mxy +
                       A_rho[8] * myy;

    real_t inv_rho_denom = real_t(1) / rho_denom;

    rho = rho_I * inv_rho_denom;

    const real_t ur = ux;
    const real_t ut = uy;
    const real_t mrr = mxx;
    const real_t mrt = mxy;
    const real_t mtt = myy;

    ux = ur * c - ut * s;
    uy = ut * c + ur * s;

    mxx = c * c * mrr - r::two * c * s * mrt + s * s * mtt;
    mxy = c * s * mrr + (c * c - s * s) * mrt - c * s * mtt;
    myy = s * s * mrr + r::two * c * s * mrt + c * c * mtt;
}