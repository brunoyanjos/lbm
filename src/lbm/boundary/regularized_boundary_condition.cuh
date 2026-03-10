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

__device__ __forceinline__ void apply_boundary_dirichlet(
    real_t *__restrict__ pop,
    uint32_t valid_mask,
    real_t &rho,
    real_t ux, real_t uy,
    real_t &mxx, real_t &mxy, real_t &myy)
{
    const uint32_t outgoing_mask = valid_mask;
    const uint32_t incoming_mask = mask_opp(valid_mask);

    // soma das válidas
    real_t rho_I = r::zero;

    real_t mxx_I = r::zero;
    real_t mxy_I = r::zero;
    real_t myy_I = r::zero;

    // from rhoI equation
    real_t sum_Os_wi = r::zero;

    real_t sum_Os_wi_Hx = r::zero;
    real_t sum_Os_wi_Hy = r::zero;

    real_t sum_Os_wi_Hxx = r::zero;
    real_t sum_Os_wi_Hxy = r::zero;
    real_t sum_Os_wi_Hyy = r::zero;

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
    A(0, 0) = sum_Is_wi_Hxx_Hxx - (r::one - Geometry::OMEGA) * sum_Os_wi_Hxx * mxx_I;
    A(0, 1) = r::two * (sum_Is_wi_Hxy_Hxx - (r::one - Geometry::OMEGA) * sum_Os_wi_Hxy * mxx_I);
    A(0, 2) = sum_Is_wi_Hyy_Hxx - (r::one - Geometry::OMEGA) * sum_Os_wi_Hyy * mxx_I;
    b[0] = sum_Os_wi * mxx_I - sum_Is_wi_Hxx +
           ux * (sum_Os_wi_Hx * mxx_I - sum_Is_wi_Hx_Hxx) + uy * (sum_Os_wi_Hy * mxx_I - sum_Is_wi_Hy_Hxx) +
           Geometry::OMEGA * ux * ux * sum_Os_wi_Hxx * mxx_I +
           r::two * Geometry::OMEGA * ux * uy * sum_Os_wi_Hxy * mxx_I +
           Geometry::OMEGA * uy * uy * sum_Os_wi_Hyy * mxx_I;

    // mxy
    A(1, 0) = sum_Is_wi_Hxx_Hxy - (r::one - Geometry::OMEGA) * sum_Os_wi_Hxx * mxy_I;
    A(1, 1) = r::two * (sum_Is_wi_Hxy_Hxy - (r::one - Geometry::OMEGA) * sum_Os_wi_Hxy * mxy_I);
    A(1, 2) = sum_Is_wi_Hyy_Hxy - (r::one - Geometry::OMEGA) * sum_Os_wi_Hyy * mxy_I;
    b[1] = sum_Os_wi * mxy_I - sum_Is_wi_Hxy +
           ux * (sum_Os_wi_Hx * mxy_I - sum_Is_wi_Hx_Hxy) + uy * (sum_Os_wi_Hy * mxy_I - sum_Is_wi_Hy_Hxy) +
           Geometry::OMEGA * ux * ux * sum_Os_wi_Hxx * mxy_I +
           r::two * Geometry::OMEGA * ux * uy * sum_Os_wi_Hxy * mxy_I +
           Geometry::OMEGA * uy * uy * sum_Os_wi_Hyy * mxy_I;

    // myy
    A(2, 0) = sum_Is_wi_Hxx_Hyy - (r::one - Geometry::OMEGA) * sum_Os_wi_Hxx * myy_I;
    A(2, 1) = r::two * (sum_Is_wi_Hxy_Hyy - (r::one - Geometry::OMEGA) * sum_Os_wi_Hxy * myy_I);
    A(2, 2) = sum_Is_wi_Hyy_Hyy - (r::one - Geometry::OMEGA) * sum_Os_wi_Hyy * myy_I;
    b[2] = sum_Os_wi * myy_I - sum_Is_wi_Hyy +
           ux * (sum_Os_wi_Hx * myy_I - sum_Is_wi_Hx_Hyy) + uy * (sum_Os_wi_Hy * myy_I - sum_Is_wi_Hy_Hyy) +
           Geometry::OMEGA * ux * ux * sum_Os_wi_Hxx * myy_I +
           r::two * Geometry::OMEGA * ux * uy * sum_Os_wi_Hyy * myy_I +
           Geometry::OMEGA * uy * uy * sum_Os_wi_Hyy * myy_I;

    gaussianElimination<3>(Aflat, b, m);

    mxx = m[0];
    mxy = m[1];
    myy = m[2];

    const real_t rho_denominator = sum_Os_wi +
                                   sum_Os_wi_Hx * ux + sum_Os_wi_Hy * uy +
                                   (r::one - Geometry::OMEGA) * sum_Os_wi_Hxx * mxx + r::two * (r::one - Geometry::OMEGA) * sum_Os_wi_Hxy * mxy + (r::one - Geometry::OMEGA) * sum_Os_wi_Hyy * myy +
                                   Geometry::OMEGA * sum_Os_wi_Hxx * ux * ux + r::two * Geometry::OMEGA * sum_Os_wi_Hxy * ux * uy + Geometry::OMEGA * sum_Os_wi_Hyy * uy * uy;
    const real_t inv_rho = real_t(1) / rho_denominator;

    rho = rho_I * inv_rho;
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