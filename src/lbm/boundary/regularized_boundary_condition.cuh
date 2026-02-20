#pragma once
#include "../../core/types.cuh"
#include "../../core/physics.h"
#include "../../core/math_utils.cuh"
#include "../../core/linear_solver.cuh"
#include "../stencil_active.cuh"
#include "../domain/mask_utils.cuh"
#include "../state/lbm_state.cuh"
#include <cstdint>
#include <cstdio>

__device__ __forceinline__ void apply_boundary(
    real_t *__restrict__ pop,
    uint32_t valid_mask,
    real_t &rho,
    real_t ux_bc, real_t uy_bc,
    real_t &mxx, real_t &mxy, real_t &myy, int x, int y, uint8_t wall_id)
{
    const uint32_t outgoing_mask = valid_mask;
    const uint32_t incoming_mask = mask_opp(valid_mask);

    // soma das válidas
    real_t rho_I = real_t(0);

    real_t mxx_I = real_t(0);
    real_t mxy_I = real_t(0);
    real_t myy_I = real_t(0);

    real_t A = real_t(0);

    real_t Bxx = real_t(0);
    real_t Bxy = real_t(0);
    real_t Byy = real_t(0);

    real_t A_Hxx = real_t(0);
    real_t A_Hxy = real_t(0);
    real_t A_Hyy = real_t(0);

    real_t Bxx_Hxx = real_t(0);
    real_t Bxx_Hxy = real_t(0);
    real_t Bxx_Hyy = real_t(0);

    real_t Bxy_Hxx = real_t(0);
    real_t Bxy_Hxy = real_t(0);
    real_t Bxy_Hyy = real_t(0);

    real_t Byy_Hxx = real_t(0);
    real_t Byy_Hxy = real_t(0);
    real_t Byy_Hyy = real_t(0);

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        real_t cx, cy, Hxx, Hxy, Hyy;
        Stencil::basis2(i, cx, cy, Hxx, Hxy, Hyy);

        const real_t A_i = Stencil::w(i) * (real_t(1) + Stencil::as2 * ux_bc * cx + Stencil::as2 * uy_bc * cy);
        const real_t Bxx_i = Stencil::w(i) * Stencil::as4 * real_t(0.5) * Hxx;
        const real_t Bxy_i = Stencil::w(i) * Stencil::as4 * real_t(0.5) * Hxy;
        const real_t Byy_i = Stencil::w(i) * Stencil::as4 * real_t(0.5) * Hyy;

        if (dir_valid(incoming_mask, i))
        {
            rho_I += pop[i];

            mxx_I += pop[i] * Hxx;
            mxy_I += pop[i] * Hxy;
            myy_I += pop[i] * Hyy;

            A_Hxx += A_i * Hxx;
            A_Hxy += A_i * Hxy;
            A_Hyy += A_i * Hyy;

            Bxx_Hxx += Bxx_i * Hxx;
            Bxy_Hxx += Bxy_i * Hxx;
            Byy_Hxx += Byy_i * Hxx;

            Bxx_Hxy += Bxx_i * Hxy;
            Bxy_Hxy += Bxy_i * Hxy;
            Byy_Hxy += Byy_i * Hxy;

            Bxx_Hyy += Bxx_i * Hyy;
            Bxy_Hyy += Bxy_i * Hyy;
            Byy_Hyy += Byy_i * Hyy;
        }

        if (dir_valid(outgoing_mask, i))
        {
            A += A_i;

            Bxx += Bxx_i;
            Bxy += Bxy_i;
            Byy += Byy_i;
        }
    }

    const real_t inv_rho_I = real_t(1) / rho_I;

    mxx_I *= inv_rho_I;
    mxy_I *= inv_rho_I;
    myy_I *= inv_rho_I;

    const real_t u_sum = ux_bc * ux_bc * Bxx +
                         real_t(2) * ux_bc * uy_bc * Bxy +
                         uy_bc * uy_bc * Byy;

    constexpr int Nsys = 3;

    real_t Uflat[Nsys * Nsys] = {
        real_t(0), real_t(0), real_t(0),
        real_t(0), real_t(0), real_t(0),
        real_t(0), real_t(0), real_t(0)};

    real_t d[Nsys] = {real_t(0), real_t(0), real_t(0)};
    real_t m[Nsys] = {real_t(0), real_t(0), real_t(0)};

    auto U = [&](int r, int c) -> real_t &
    { return Uflat[r * Nsys + c]; };

    // mxx
    U(0, 0) = (real_t(1) - OMEGA) * Bxx * mxx_I - Bxx_Hxx;
    U(0, 1) = real_t(2) * ((real_t(1) - OMEGA) * Bxy * mxx_I - Bxy_Hxx);
    U(0, 2) = (real_t(1) - OMEGA) * Byy * mxx_I - Byy_Hxx;
    d[0] = A_Hxx - (A + OMEGA * u_sum) * mxx_I;

    // mxy
    U(1, 0) = (real_t(1) - OMEGA) * Bxx * mxy_I - Bxx_Hxy;
    U(1, 1) = real_t(2) * ((real_t(1) - OMEGA) * Bxy * mxy_I - Bxy_Hxy);
    U(1, 2) = (real_t(1) - OMEGA) * Byy * mxy_I - Byy_Hxy;
    d[1] = A_Hxy - (A + OMEGA * u_sum) * mxy_I;

    // myy
    U(2, 0) = (real_t(1) - OMEGA) * Bxx * myy_I - Bxx_Hyy;
    U(2, 1) = real_t(2) * ((real_t(1) - OMEGA) * Bxy * myy_I - Bxy_Hyy);
    U(2, 2) = (real_t(1) - OMEGA) * Byy * myy_I - Byy_Hyy;
    d[2] = A_Hyy - (A + OMEGA * u_sum) * myy_I;

    gaussianElimination<3>(Uflat, d, m);

    mxx = m[0];
    mxy = m[1];
    myy = m[2];

    // Uflat é row-major: U(r,c) = Uflat[r*3+c]
    // const real_t a11 = Uflat[0], a12 = Uflat[1], a13 = Uflat[2];
    // const real_t a21 = Uflat[3], a22 = Uflat[4], a23 = Uflat[5];
    // const real_t a31 = Uflat[6], a32 = Uflat[7], a33 = Uflat[8];

    // const real_t b1 = d[0], b2 = d[1], b3 = d[2];

    // const real_t den =
    //     a13 * a22 * a31 - a12 * a23 * a31 - a13 * a21 * a32 + a11 * a23 * a32 + a12 * a21 * a33 - a11 * a22 * a33;

    // const real_t inv_den = real_t(1) / den;

    // // exatamente como no código que você mandou (atenção ao sinal de mxy)
    // mxx = (a23 * a32 * b1 - a22 * a33 * b1 - a13 * a32 * b2 + a12 * a33 * b2 + a13 * a22 * b3 - a12 * a23 * b3) * inv_den;
    // mxy = -(a23 * a31 * b1 - a21 * a33 * b1 - a13 * a31 * b2 + a11 * a33 * b2 + a13 * a21 * b3 - a11 * a23 * b3) * inv_den;
    // myy = (a22 * a31 * b1 - a21 * a32 * b1 - a12 * a31 * b2 + a11 * a32 * b2 + a12 * a21 * b3 - a11 * a22 * b3) * inv_den;

    const real_t mom_sum = mxx * Bxx + real_t(2) * mxy * Bxy + myy * Byy;

    const real_t rho_denominator = A + (real_t(1) - OMEGA) * mom_sum + OMEGA * u_sum;
    const real_t inv_rho = real_t(1) / rho_denominator;

    rho = rho_I * inv_rho;
}
