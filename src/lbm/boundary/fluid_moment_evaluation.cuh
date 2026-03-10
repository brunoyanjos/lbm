#pragma once

#include "../../core/types.cuh"
#include "../../core/linear_solver.cuh"
#include "../domain/mask_utils.cuh"
#include "../../geometries/active_geometry.cuh"
#include <cstdint>

__device__ inline void evaluate_fluid_node(
    real_t *__restrict__ pop,
    uint32_t valid_mask,
    real_t &rho,
    real_t &ux, real_t &uy,
    real_t &mxx, real_t &mxy, real_t &myy)
{
    const uint32_t outgoing_mask = valid_mask;
    const uint32_t incoming_mask = mask_opp(valid_mask);

    real_t A_coeff[40]{};
    real_t b_coeff[5]{};

    // incoming moments
    real_t rho_I = real_t(0);

    real_t ux_I = real_t(0);
    real_t uy_I = real_t(0);

    real_t mxx_I = real_t(0);
    real_t mxy_I = real_t(0);
    real_t myy_I = real_t(0);

    // from rho equation
    real_t sum_Os_wi = real_t(0);

    real_t sum_Os_wi_Hx = real_t(0);
    real_t sum_Os_wi_Hy = real_t(0);

    real_t sum_Os_wi_Hxx = real_t(0);
    real_t sum_Os_wi_Hxy = real_t(0);
    real_t sum_Os_wi_Hyy = real_t(0);

    // from uxI equation
    real_t sum_Is_wi_Hx = real_t(0);

    real_t sum_Is_wi_Hx_Hx = real_t(0);
    real_t sum_Is_wi_Hy_Hx = real_t(0);

    real_t sum_Is_wi_Hxx_Hx = real_t(0);
    real_t sum_Is_wi_Hxy_Hx = real_t(0);
    real_t sum_Is_wi_Hyy_Hx = real_t(0);

    // from uyI equation
    real_t sum_Is_wi_Hy = real_t(0);

    real_t sum_Is_wi_Hx_Hy = real_t(0);
    real_t sum_Is_wi_Hy_Hy = real_t(0);

    real_t sum_Is_wi_Hxx_Hy = real_t(0);
    real_t sum_Is_wi_Hxy_Hy = real_t(0);
    real_t sum_Is_wi_Hyy_Hy = real_t(0);

    // from mxxI equation
    real_t sum_Is_wi_Hxx = real_t(0);

    real_t sum_Is_wi_Hx_Hxx = real_t(0);
    real_t sum_Is_wi_Hy_Hxx = real_t(0);

    real_t sum_Is_wi_Hxx_Hxx = real_t(0);
    real_t sum_Is_wi_Hxy_Hxx = real_t(0);
    real_t sum_Is_wi_Hyy_Hxx = real_t(0);

    // from mxyI equation
    real_t sum_Is_wi_Hxy = real_t(0);

    real_t sum_Is_wi_Hx_Hxy = real_t(0);
    real_t sum_Is_wi_Hy_Hxy = real_t(0);

    real_t sum_Is_wi_Hxx_Hxy = real_t(0);
    real_t sum_Is_wi_Hxy_Hxy = real_t(0);
    real_t sum_Is_wi_Hyy_Hxy = real_t(0);

    // from myyI equation
    real_t sum_Is_wi_Hyy = real_t(0);

    real_t sum_Is_wi_Hx_Hyy = real_t(0);
    real_t sum_Is_wi_Hy_Hyy = real_t(0);

    real_t sum_Is_wi_Hxx_Hyy = real_t(0);
    real_t sum_Is_wi_Hxy_Hyy = real_t(0);
    real_t sum_Is_wi_Hyy_Hyy = real_t(0);

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        real_t cx, cy, Hxx, Hxy, Hyy;
        Stencil::basis2(i, cx, cy, Hxx, Hxy, Hyy);

        const real_t Hx_factor = Stencil::w(i) * Stencil::as2 * cx;
        const real_t Hy_factor = Stencil::w(i) * Stencil::as2 * cy;

        const real_t Hxx_factor = Stencil::w(i) * Stencil::as4 * real_t(0.5) * Hxx;
        const real_t Hxy_factor = Stencil::w(i) * Stencil::as4 * real_t(0.5) * Hxy;
        const real_t Hyy_factor = Stencil::w(i) * Stencil::as4 * real_t(0.5) * Hyy;

        if (dir_valid(incoming_mask, i))
        {
            rho_I += pop[i];
            ux_I += pop[i] * cx;
            uy_I += pop[i] * cy;
            mxx_I += pop[i] * Hxx;
            mxy_I += pop[i] * Hxy;
            myy_I += pop[i] * Hyy;

            sum_Is_wi_Hx += Stencil::w(i) * cx;
            sum_Is_wi_Hx_Hx += Hx_factor * cx;
            sum_Is_wi_Hy_Hx += Hy_factor * cx;
            sum_Is_wi_Hxx_Hx += Hxx_factor * cx;
            sum_Is_wi_Hxy_Hx += Hxy_factor * cx;
            sum_Is_wi_Hyy_Hx += Hyy_factor * cx;

            sum_Is_wi_Hy += Stencil::w(i) * cy;
            sum_Is_wi_Hx_Hy += Hx_factor * cy;
            sum_Is_wi_Hy_Hy += Hy_factor * cy;
            sum_Is_wi_Hxx_Hy += Hxx_factor * cy;
            sum_Is_wi_Hxy_Hy += Hxy_factor * cy;
            sum_Is_wi_Hyy_Hy += Hyy_factor * cy;

            sum_Is_wi_Hxx += Stencil::w(i) * Hxx;
            sum_Is_wi_Hx_Hxx += Hx_factor * Hxx;
            sum_Is_wi_Hy_Hxx += Hy_factor * Hxx;
            sum_Is_wi_Hxx_Hxx += Hxx_factor * Hxx;
            sum_Is_wi_Hxy_Hxx += Hxy_factor * Hxx;
            sum_Is_wi_Hyy_Hxx += Hyy_factor * Hxx;

            sum_Is_wi_Hxy += Stencil::w(i) * Hxy;
            sum_Is_wi_Hx_Hxy += Hx_factor * Hxy;
            sum_Is_wi_Hy_Hxy += Hy_factor * Hxy;
            sum_Is_wi_Hxx_Hxy += Hxx_factor * Hxy;
            sum_Is_wi_Hxy_Hxy += Hxy_factor * Hxy;
            sum_Is_wi_Hyy_Hxy += Hyy_factor * Hxy;

            sum_Is_wi_Hyy += Stencil::w(i) * Hyy;
            sum_Is_wi_Hx_Hyy += Hx_factor * Hyy;
            sum_Is_wi_Hy_Hyy += Hy_factor * Hyy;
            sum_Is_wi_Hxx_Hyy += Hxx_factor * Hyy;
            sum_Is_wi_Hxy_Hyy += Hxy_factor * Hyy;
            sum_Is_wi_Hyy_Hyy += Hyy_factor * Hyy;
        }

        if (dir_valid(outgoing_mask, i))
        {
            sum_Os_wi += Stencil::w(i);
            sum_Os_wi_Hx += Hx_factor;
            sum_Os_wi_Hy += Hy_factor;
            sum_Os_wi_Hxx += Hxx_factor;
            sum_Os_wi_Hxy += Hxy_factor;
            sum_Os_wi_Hyy += Hyy_factor;
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
    A(0, 0) = sum_Is_wi_Hx_Hx - sum_Os_wi_Hx * ux_I;                                                 // ux
    A(0, 1) = sum_Is_wi_Hy_Hx - sum_Os_wi_Hy * ux_I;                                                 // uy
    A(0, 2) = -Geometry::OMEGA * sum_Os_wi_Hxx * ux_I;                                               // ux ux
    A(0, 3) = -real_t(2) * Geometry::OMEGA * sum_Os_wi_Hxy * ux_I;                                   // ux uy
    A(0, 4) = -Geometry::OMEGA * sum_Os_wi_Hyy * ux_I;                                               // uy uy
    A(0, 5) = sum_Is_wi_Hxx_Hx - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hxx * ux_I;               // mxx
    A(0, 6) = real_t(2) * (sum_Is_wi_Hxy_Hx - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hxy * ux_I); // mxy
    A(0, 7) = sum_Is_wi_Hyy_Hx - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hyy * ux_I;               // myy

    b_coeff[0] = ux_I * sum_Os_wi - sum_Is_wi_Hx;

    // uyI equation
    A(1, 0) = sum_Is_wi_Hx_Hy - sum_Os_wi_Hx * uy_I;                                                 // ux
    A(1, 1) = sum_Is_wi_Hy_Hy - sum_Os_wi_Hy * uy_I;                                                 // uy
    A(1, 2) = -Geometry::OMEGA * sum_Os_wi_Hxx * uy_I;                                               // ux ux
    A(1, 3) = -real_t(2) * Geometry::OMEGA * sum_Os_wi_Hxy * uy_I;                                   // ux uy
    A(1, 4) = -Geometry::OMEGA * sum_Os_wi_Hyy * uy_I;                                               // uy uy
    A(1, 5) = sum_Is_wi_Hxx_Hy - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hxx * uy_I;               // mxx
    A(1, 6) = real_t(2) * (sum_Is_wi_Hxy_Hy - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hxy * uy_I); // mxy
    A(1, 7) = sum_Is_wi_Hyy_Hy - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hyy * uy_I;               // myy

    b_coeff[1] = uy_I * sum_Os_wi - sum_Is_wi_Hy;

    // mxxI equation
    A(2, 0) = sum_Is_wi_Hx_Hxx - sum_Os_wi_Hx * mxx_I;                                                 // ux
    A(2, 1) = sum_Is_wi_Hy_Hxx - sum_Os_wi_Hy * mxx_I;                                                 // uy
    A(2, 2) = -Geometry::OMEGA * sum_Os_wi_Hxx * mxx_I;                                                // ux ux
    A(2, 3) = -real_t(2) * Geometry::OMEGA * sum_Os_wi_Hxy * mxx_I;                                    // ux uy
    A(2, 4) = -Geometry::OMEGA * sum_Os_wi_Hyy * mxx_I;                                                // uy uy
    A(2, 5) = sum_Is_wi_Hxx_Hxx - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hxx * mxx_I;               // mxx
    A(2, 6) = real_t(2) * (sum_Is_wi_Hxy_Hxx - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hxy * mxx_I); // mxy
    A(2, 7) = sum_Is_wi_Hyy_Hxx - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hyy * mxx_I;               // myy

    b_coeff[2] = mxx_I * sum_Os_wi - sum_Is_wi_Hxx;

    // mxyI equation
    A(3, 0) = sum_Is_wi_Hx_Hxy - sum_Os_wi_Hx * mxy_I;                                                 // ux
    A(3, 1) = sum_Is_wi_Hy_Hxy - sum_Os_wi_Hy * mxy_I;                                                 // uy
    A(3, 2) = -Geometry::OMEGA * sum_Os_wi_Hxx * mxy_I;                                                // ux ux
    A(3, 3) = -real_t(2) * Geometry::OMEGA * sum_Os_wi_Hxy * mxy_I;                                    // ux uy
    A(3, 4) = -Geometry::OMEGA * sum_Os_wi_Hyy * mxy_I;                                                // uy uy
    A(3, 5) = sum_Is_wi_Hxx_Hxy - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hxx * mxy_I;               // mxx
    A(3, 6) = real_t(2) * (sum_Is_wi_Hxy_Hxy - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hxy * mxy_I); // mxy
    A(3, 7) = sum_Is_wi_Hyy_Hxy - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hyy * mxy_I;               // myy

    b_coeff[3] = mxy_I * sum_Os_wi - sum_Is_wi_Hxy;

    // myyI equation
    A(4, 0) = sum_Is_wi_Hx_Hyy - sum_Os_wi_Hx * myy_I;                                                 // ux
    A(4, 1) = sum_Is_wi_Hy_Hyy - sum_Os_wi_Hy * myy_I;                                                 // uy
    A(4, 2) = -Geometry::OMEGA * sum_Os_wi_Hxx * myy_I;                                                // ux ux
    A(4, 3) = -real_t(2) * Geometry::OMEGA * sum_Os_wi_Hxy * myy_I;                                    // ux uy
    A(4, 4) = -Geometry::OMEGA * sum_Os_wi_Hyy * myy_I;                                                // uy uy
    A(4, 5) = sum_Is_wi_Hxx_Hyy - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hxx * myy_I;               // mxx
    A(4, 6) = real_t(2) * (sum_Is_wi_Hxy_Hyy - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hxy * myy_I); // mxy
    A(4, 7) = sum_Is_wi_Hyy_Hyy - (real_t(1) - Geometry::OMEGA) * sum_Os_wi_Hyy * myy_I;               // myy

    b_coeff[4] = myy_I * sum_Os_wi - sum_Is_wi_Hyy;

    real_t A_gauss[5 * 5]{};
    real_t b_gauss[5]{};
    real_t moms[5]{};

    auto U = [&](int r, int c) -> real_t &
    { return A_gauss[r * 5 + c]; };

    real_t error = real_t(1);
    int it = 0;
    const int it_max = 50;

    while (error > real_t(1e-6) && it++ < it_max)
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

    real_t rho_denom = (sum_Os_wi +
                        sum_Os_wi_Hx * ux +
                        sum_Os_wi_Hy * uy +
                        sum_Os_wi_Hxx * ((1 - Geometry::OMEGA) * mxx + Geometry::OMEGA * ux * ux) +
                        real_t(2) * sum_Os_wi_Hxy * ((1 - Geometry::OMEGA) * mxy + Geometry::OMEGA * ux * uy) +
                        sum_Os_wi_Hyy * ((1 - Geometry::OMEGA) * myy + Geometry::OMEGA * uy * uy));

    real_t inv_rho_denom = real_t(1) / rho_denom;

    rho = rho_I * inv_rho_denom;
}