#pragma once
#include "../../core/types.cuh"
#include "../../geometries/active_geometry.cuh"
#include "../../core/math_utils.cuh"
#include "../../core/linear_solver.cuh"
#include "../stencil_active.cuh"
#include "../domain/mask_utils.cuh"
#include "../state/lbm_state.cuh"
#include <cstdint>
#include <cstdio>

struct SystemVariable
{
    real_t mom_I[4]{};

    real_t rho[9]{};

    real_t mxx[6]{};
    real_t mxy[6]{};
    real_t myy[6]{};

    __device__ __forceinline__ void normalize_moments()
    {
        const real_t inv_rho = r::one / mom_I[0];

        for (int i = 1; i < 4; ++i)
        {
            mom_I[i] *= inv_rho;
        }
    }
};

__device__ __forceinline__ void
apply_boundary_dirichlet(
    real_t *__restrict__ pop,
    uint32_t valid_mask,
    real_t &rho,
    real_t ux, real_t uy,
    real_t &mxx, real_t &mxy, real_t &myy)
{
    const uint32_t outgoing_mask = valid_mask;
    const uint32_t incoming_mask = mask_opp(valid_mask);

    SystemVariable Am{};

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        real_t cx, cy, Hxx, Hxy, Hyy;
        Stencil::basis2(i, cx, cy, Hxx, Hxy, Hyy);

        const real_t factor[6] = {
            Stencil::w(i),
            Stencil::w(i) * Stencil::as2 * cx,
            Stencil::w(i) * Stencil::as2 * cy,
            Stencil::w(i) * Stencil::as4 * r::half * Hxx,
            Stencil::w(i) * Stencil::as4 * Hxy,
            Stencil::w(i) * Stencil::as4 * r::half * Hyy};

        if (dir_valid(incoming_mask, i))
        {
            Am.mom_I[0] += pop[i];

            Am.mom_I[1] += pop[i] * Hxx;
            Am.mom_I[2] += pop[i] * Hxy;
            Am.mom_I[3] += pop[i] * Hyy;

#pragma unroll
            for (int j = 0; j < 6; ++j)
            {
                Am.mxx[j] += factor[j] * Hxx;
                Am.mxy[j] += factor[j] * Hxy;
                Am.myy[j] += factor[j] * Hyy;
            }
        }

        if (dir_valid(outgoing_mask, i))
        {
            Am.rho[0] += factor[0];
            Am.rho[1] += factor[1];
            Am.rho[2] += factor[2];
            Am.rho[3] += (r::one - Geometry::OMEGA) * factor[3];
            Am.rho[4] += (r::one - Geometry::OMEGA) * factor[4];
            Am.rho[5] += (r::one - Geometry::OMEGA) * factor[5];
            Am.rho[6] += Geometry::OMEGA * factor[3];
            Am.rho[7] += Geometry::OMEGA * factor[4];
            Am.rho[8] += Geometry::OMEGA * factor[5];
        }
    }

    Am.normalize_moments();

    constexpr int Nsys = 3;

    real_t Aflat[Nsys * Nsys]{};

    real_t b[Nsys]{};
    real_t m[Nsys]{};

    auto A = [&](int r, int c) -> real_t &
    { return Aflat[r * Nsys + c]; };

    A(0, 0) = Am.mxx[3] - Am.rho[3] * Am.mom_I[1];
    A(0, 1) = Am.mxx[4] - Am.rho[4] * Am.mom_I[1];
    A(0, 2) = Am.mxx[5] - Am.rho[5] * Am.mom_I[1];
    b[0] = Am.mom_I[1] * (Am.rho[0] +
                          ux * Am.rho[1] + uy * Am.rho[2] +
                          ux * ux * Am.rho[6] + ux * uy * Am.rho[7] + uy * uy * Am.rho[8]) -
           (Am.mxx[0] + ux * Am.mxx[1] + uy * Am.mxx[2]);

    A(1, 0) = Am.mxy[3] - Am.rho[3] * Am.mom_I[2];
    A(1, 1) = Am.mxy[4] - Am.rho[4] * Am.mom_I[2];
    A(1, 2) = Am.mxy[5] - Am.rho[5] * Am.mom_I[2];
    b[1] = Am.mom_I[2] * (Am.rho[0] +
                          ux * Am.rho[1] + uy * Am.rho[2] +
                          ux * ux * Am.rho[6] + ux * uy * Am.rho[7] + uy * uy * Am.rho[8]) -
           (Am.mxy[0] + ux * Am.mxy[1] + uy * Am.mxy[2]);

    A(2, 0) = Am.myy[3] - Am.rho[3] * Am.mom_I[3];
    A(2, 1) = Am.myy[4] - Am.rho[4] * Am.mom_I[3];
    A(2, 2) = Am.myy[5] - Am.rho[5] * Am.mom_I[3];
    b[2] = Am.mom_I[3] * (Am.rho[0] +
                          ux * Am.rho[1] + uy * Am.rho[2] +
                          ux * ux * Am.rho[6] + ux * uy * Am.rho[7] + uy * uy * Am.rho[8]) -
           (Am.myy[0] + ux * Am.myy[1] + uy * Am.myy[2]);

    gaussianElimination<3>(Aflat, b, m);

    mxx = m[0];
    mxy = m[1];
    myy = m[2];

    const real_t rho_denominator = Am.rho[0] + ux * Am.rho[1] + uy * Am.rho[2] +
                                   mxx * Am.rho[3] + mxy * Am.rho[4] + myy * Am.rho[5] +
                                   ux * ux * Am.rho[6] + ux * uy * Am.rho[7] + uy * uy * Am.rho[8];
    const real_t inv_rho = r::one / rho_denominator;

    rho = Am.mom_I[0] * inv_rho;
}

__device__ __forceinline__ void apply_boundary_inlet(
    real_t *__restrict__ pop,
    uint32_t valid_mask,
    real_t &rho,
    real_t ux, real_t uy,
    real_t &mxx, real_t &mxy, real_t &myy)
{
    const uint32_t incoming_mask = mask_opp(valid_mask);

    // soma das válidas
    real_t rho_I = r::zero;

    real_t mxx_I = r::zero;
    real_t mxy_I = r::zero;
    real_t myy_I = r::zero;

    // from rhoI equation
    real_t sum_Is_wi = r::zero;

    real_t sum_Is_wi_Hx = r::zero;
    real_t sum_Is_wi_Hy = r::zero;

    // from mxxI equation
    real_t sum_Is_wi_Hxx = r::zero;

    real_t sum_Is_wi_Hx_Hxx = r::zero;
    real_t sum_Is_wi_Hy_Hxx = r::zero;

    real_t sum_Is_wi_Hxx_Hxx = r::zero;
    real_t sum_Is_wi_Hxy_Hxx = r::zero;
    real_t sum_Is_wi_Hyy_Hxx = r::zero;

    // from mxyI equation
    real_t sum_Is_wi_Hxy = r::zero;

    real_t sum_Is_wi_Hx_Hxy = r::zero;
    real_t sum_Is_wi_Hy_Hxy = r::zero;

    real_t sum_Is_wi_Hxx_Hxy = r::zero;
    real_t sum_Is_wi_Hxy_Hxy = r::zero;
    real_t sum_Is_wi_Hyy_Hxy = r::zero;

    // from myyI equation
    real_t sum_Is_wi_Hyy = r::zero;

    real_t sum_Is_wi_Hx_Hyy = r::zero;
    real_t sum_Is_wi_Hy_Hyy = r::zero;

    real_t sum_Is_wi_Hxx_Hyy = r::zero;
    real_t sum_Is_wi_Hxy_Hyy = r::zero;
    real_t sum_Is_wi_Hyy_Hyy = r::zero;

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        real_t cx, cy, Hxx, Hxy, Hyy;
        Stencil::basis2(i, cx, cy, Hxx, Hxy, Hyy);

        const real_t Hx_factor = Stencil::w(i) * Stencil::as2 * cx;
        const real_t Hy_factor = Stencil::w(i) * Stencil::as2 * cy;

        const real_t Hxx_factor = Stencil::w(i) * Stencil::as4 * r::half * Hxx;
        const real_t Hxy_factor = Stencil::w(i) * Stencil::as4 * r::half * Hxy;
        const real_t Hyy_factor = Stencil::w(i) * Stencil::as4 * r::half * Hyy;

        if (dir_valid(incoming_mask, i))
        {
            rho_I += pop[i];

            mxx_I += pop[i] * Hxx;
            mxy_I += pop[i] * Hxy;
            myy_I += pop[i] * Hyy;

            sum_Is_wi += Stencil::w(i);

            sum_Is_wi_Hx += Hx_factor;
            sum_Is_wi_Hy += Hy_factor;

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
    }

    const real_t inv_rho_I = real_t(1) / rho_I;

    mxx_I *= inv_rho_I;
    mxy_I *= inv_rho_I;
    myy_I *= inv_rho_I;

    constexpr int Nsys = 3;

    real_t Aflat[Nsys * Nsys]{};

    real_t b[Nsys]{};
    real_t m[Nsys]{};

    auto A = [&](int r, int c) -> real_t &
    { return Aflat[r * Nsys + c]; };

    // mxx
    A(0, 0) = sum_Is_wi_Hxx_Hxx - Stencil::as4 * r::half * sum_Is_wi_Hxx * mxx_I;
    A(0, 1) = r::two * (sum_Is_wi_Hxy_Hxx - Stencil::as4 * r::half * sum_Is_wi_Hxy * mxx_I);
    A(0, 2) = sum_Is_wi_Hyy_Hxx - Stencil::as4 * r::half * sum_Is_wi_Hyy * mxx_I;
    b[0] = sum_Is_wi * mxx_I - sum_Is_wi_Hxx +
           ux * (sum_Is_wi_Hx * mxx_I - sum_Is_wi_Hx_Hxx) +
           uy * (sum_Is_wi_Hy * mxx_I - sum_Is_wi_Hy_Hxx);

    // mxy
    A(1, 0) = sum_Is_wi_Hxx_Hxx - Stencil::as4 * r::half * sum_Is_wi_Hxx * mxx_I;
    A(1, 1) = r::two * (sum_Is_wi_Hxy_Hxx - Stencil::as4 * r::half * sum_Is_wi_Hxy * mxx_I);
    A(1, 2) = sum_Is_wi_Hyy_Hxx - Stencil::as4 * r::half * sum_Is_wi_Hyy * mxx_I;
    b[1] = sum_Is_wi * mxy_I - sum_Is_wi_Hxy +
           ux * (sum_Is_wi_Hx * mxy_I - sum_Is_wi_Hx_Hxy) +
           uy * (sum_Is_wi_Hy * mxy_I - sum_Is_wi_Hy_Hxy);

    // myy
    A(2, 0) = sum_Is_wi_Hxx_Hxx - Stencil::as4 * r::half * sum_Is_wi_Hxx * mxx_I;
    A(2, 1) = r::two * (sum_Is_wi_Hxy_Hxx - Stencil::as4 * r::half * sum_Is_wi_Hxy * mxx_I);
    A(2, 2) = sum_Is_wi_Hyy_Hxx - Stencil::as4 * r::half * sum_Is_wi_Hyy * mxx_I;
    b[2] = sum_Is_wi * myy_I - sum_Is_wi_Hyy +
           ux * (sum_Is_wi_Hx * myy_I - sum_Is_wi_Hx_Hyy) +
           uy * (sum_Is_wi_Hy * myy_I - sum_Is_wi_Hy_Hyy);

    gaussianElimination<3>(Aflat, b, m);

    mxx = m[0];
    mxy = m[1];
    myy = m[2];

    const real_t rho_denominator = sum_Is_wi +
                                   sum_Is_wi_Hx * ux + sum_Is_wi_Hy * uy +
                                   Stencil::as2 * r::half * sum_Is_wi_Hxx * mxx +
                                   Stencil::as2 * sum_Is_wi_Hxy * mxy +
                                   Stencil::as2 * r::half * sum_Is_wi_Hyy * myy;

    const real_t inv_rho = real_t(1) / rho_denominator;

    rho = rho_I * inv_rho;
}

__device__ __forceinline__ void apply_boundary_outlet(
    real_t *__restrict__ pop,
    uint32_t valid_mask,
    real_t rho,
    real_t ux, real_t uy,
    real_t &mxx, real_t &mxy, real_t &myy)
{
    const uint32_t incoming_mask = mask_opp(valid_mask);

    // soma das válidas
    real_t rho_I = r::zero;

    real_t mxx_I = r::zero;
    real_t mxy_I = r::zero;
    real_t myy_I = r::zero;

    // from mxxI equation
    real_t sum_Is_wi_Hxx = r::zero;

    real_t sum_Is_wi_Hx_Hxx = r::zero;
    real_t sum_Is_wi_Hy_Hxx = r::zero;

    real_t sum_Is_wi_Hxx_Hxx = r::zero;
    real_t sum_Is_wi_Hxy_Hxx = r::zero;
    real_t sum_Is_wi_Hyy_Hxx = r::zero;

    // from mxyI equation
    real_t sum_Is_wi_Hxy = r::zero;

    real_t sum_Is_wi_Hx_Hxy = r::zero;
    real_t sum_Is_wi_Hy_Hxy = r::zero;

    real_t sum_Is_wi_Hxx_Hxy = r::zero;
    real_t sum_Is_wi_Hxy_Hxy = r::zero;
    real_t sum_Is_wi_Hyy_Hxy = r::zero;

    // from myyI equation
    real_t sum_Is_wi_Hyy = r::zero;

    real_t sum_Is_wi_Hx_Hyy = r::zero;
    real_t sum_Is_wi_Hy_Hyy = r::zero;

    real_t sum_Is_wi_Hxx_Hyy = r::zero;
    real_t sum_Is_wi_Hxy_Hyy = r::zero;
    real_t sum_Is_wi_Hyy_Hyy = r::zero;

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        real_t cx, cy, Hxx, Hxy, Hyy;
        Stencil::basis2(i, cx, cy, Hxx, Hxy, Hyy);

        const real_t Hx_factor = Stencil::w(i) * Stencil::as2 * cx;
        const real_t Hy_factor = Stencil::w(i) * Stencil::as2 * cy;

        const real_t Hxx_factor = Stencil::w(i) * Stencil::as4 * r::half * Hxx;
        const real_t Hxy_factor = Stencil::w(i) * Stencil::as4 * r::half * Hxy;
        const real_t Hyy_factor = Stencil::w(i) * Stencil::as4 * r::half * Hyy;

        if (dir_valid(incoming_mask, i))
        {
            rho_I += pop[i];

            mxx_I += pop[i] * Hxx;
            mxy_I += pop[i] * Hxy;
            myy_I += pop[i] * Hyy;

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
    }

    const real_t inv_rho_I = real_t(1) / rho_I;

    mxx_I *= inv_rho_I;
    mxy_I *= inv_rho_I;
    myy_I *= inv_rho_I;

    constexpr int Nsys = 3;

    real_t Aflat[Nsys * Nsys]{};

    real_t b[Nsys]{};
    real_t m[Nsys]{};

    auto A = [&](int r, int c) -> real_t &
    { return Aflat[r * Nsys + c]; };

    const real_t inv_rho = r::one / rho;

    // mxx
    A(0, 0) = sum_Is_wi_Hxx_Hxx;
    A(0, 1) = r::two * sum_Is_wi_Hxy_Hxx;
    A(0, 2) = sum_Is_wi_Hyy_Hxx;
    b[0] = rho_I * inv_rho * mxx_I - sum_Is_wi_Hxx -
           ux * sum_Is_wi_Hx_Hxx - uy * sum_Is_wi_Hy_Hxx;

    // mxy
    A(1, 0) = sum_Is_wi_Hxx_Hxy;
    A(1, 1) = r::two * sum_Is_wi_Hxy_Hxy;
    A(1, 2) = sum_Is_wi_Hyy_Hxy;
    b[1] = rho_I * inv_rho * mxy_I - sum_Is_wi_Hxy -
           ux * sum_Is_wi_Hx_Hxy - uy * sum_Is_wi_Hy_Hxy;

    // myy
    A(2, 0) = sum_Is_wi_Hxx_Hyy;
    A(2, 1) = r::two * sum_Is_wi_Hxy_Hyy;
    A(2, 2) = sum_Is_wi_Hyy_Hyy;
    b[2] = rho_I * inv_rho * myy_I - sum_Is_wi_Hyy -
           ux * sum_Is_wi_Hx_Hyy - uy * sum_Is_wi_Hy_Hyy;

    gaussianElimination<3>(Aflat, b, m);

    mxx = m[0];
    mxy = m[1];
    myy = m[2];
}