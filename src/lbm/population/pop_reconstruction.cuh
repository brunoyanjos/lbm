#pragma once

#include "../../core/types.cuh"
#include "../../core/indexing.cuh"
#include "../../core/simulation_config.h"
#include "../../core/physics.h"
#include "../stencil_active.cuh"
#include "../hermite/hermite.cuh"
#include "../moment/moment_id.cuh"
#include "lbm/state/lbm_state.cuh"

template <int RegOrder, bool HighOrder>
__device__ __forceinline__ void reconstruct_streamed_pop(real_t *__restrict__ pop, real_t *__restrict__ popA, real_t *__restrict__ popB,
                                                         const LBMStateFor<RegOrder, HighOrder> &S,
                                                         int c, int x, int y)
{
#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const int cx = Stencil::cx(i);
        const int cy = Stencil::cy(i);

        const size_t n_idx = idxGlobalPeriodic(x - cx, y - cy);

        const real_t rhoA = S.d_rhoA[c][n_idx];
        const real_t rhoB = S.d_rhoB[c][n_idx];

        const real_t ux = S.d_ux[c][n_idx];
        const real_t uy = S.d_uy[c][n_idx];
        const real_t mxx = S.d_mxx[c][n_idx];
        const real_t mxy = S.d_mxy[c][n_idx];
        const real_t myy = S.d_myy[c][n_idx];

        const Normal n = S.d_n[c][n_idx];

        const real_t rho = rhoA + rhoB;
        const real_t xA = rhoA / rho;
        const real_t xB = rhoB / rho;

        const real_t MASS_DENSITY_A = MASS_A * rhoA;
        const real_t MASS_DENSITY_B = MASS_B * rhoB;

        const real_t TOTAL_MASS = MASS_DENSITY_A + MASS_DENSITY_B;

        const real_t omegaA = MASS_DENSITY_A / TOTAL_MASS;
        const real_t omegaB = MASS_DENSITY_B / TOTAL_MASS;

        const real_t uxA = ux + BETA * omegaB * n.x * Stencil::cs2;
        const real_t uyA = uy + BETA * omegaB * n.y * Stencil::cs2;

        const real_t uxB = ux - BETA * omegaA * n.x * Stencil::cs2;
        const real_t uyB = uy - BETA * omegaA * n.y * Stencil::cs2;

        const real_t THETA_A = (INV_MASS_A - r::one) * Stencil::cs2;
        const real_t THETA_B = (INV_MASS_B - r::one) * Stencil::cs2;

        const real_t fi_eq_A = Stencil::w(i) * rhoA *
                               (r::one +
                                scale_factor<MomentId::ux>() * uxA * hermite<MomentId::ux>(i) +
                                scale_factor<MomentId::uy>() * uyA * hermite<MomentId::uy>(i) +
                                scale_factor<MomentId::mxx>() * (ux * ux + THETA_A) * hermite<MomentId::mxx>(i) +
                                scale_factor<MomentId::mxy>() * ux * uy * hermite<MomentId::mxy>(i) +
                                scale_factor<MomentId::myy>() * (uy * uy + THETA_A) * hermite<MomentId::myy>(i) +
                                scale_factor<MomentId::mxxy>() * (ux * ux * uy + THETA_A * uy) * hermite<MomentId::mxxy>(i) +
                                scale_factor<MomentId::mxyy>() * (ux * uy * uy + THETA_A * ux) * hermite<MomentId::mxyy>(i) +
                                scale_factor<MomentId::mxxyy>() * (ux * ux * uy * uy + THETA_A * (ux * ux + uy * uy) + THETA_A * THETA_A) * hermite<MomentId::mxxyy>(i));

        const real_t fi_eq_B = Stencil::w(i) * rhoB *
                               (r::one +
                                scale_factor<MomentId::ux>() * uxB * hermite<MomentId::ux>(i) +
                                scale_factor<MomentId::uy>() * uyB * hermite<MomentId::uy>(i) +
                                scale_factor<MomentId::mxx>() * (ux * ux + THETA_B) * hermite<MomentId::mxx>(i) +
                                scale_factor<MomentId::mxy>() * ux * uy * hermite<MomentId::mxy>(i) +
                                scale_factor<MomentId::myy>() * (uy * uy + THETA_B) * hermite<MomentId::myy>(i) +
                                scale_factor<MomentId::mxxy>() * (ux * ux * uy + THETA_B * uy) * hermite<MomentId::mxxy>(i) +
                                scale_factor<MomentId::mxyy>() * (ux * uy * uy + THETA_B * ux) * hermite<MomentId::mxyy>(i) +
                                scale_factor<MomentId::mxxyy>() * (ux * ux * uy * uy + THETA_B * (ux * ux + uy * uy) + THETA_B * THETA_B) * hermite<MomentId::mxxyy>(i));

        const real_t fi_neq = Stencil::w(i) * TOTAL_MASS *
                              (scale_factor<MomentId::mxx>() * mxx * hermite<MomentId::mxx>(i) +
                               scale_factor<MomentId::mxy>() * mxy * hermite<MomentId::mxy>(i) +
                               scale_factor<MomentId::myy>() * myy * hermite<MomentId::myy>(i));

        popA[i] = fi_eq_A;
        popB[i] = fi_eq_B;
        pop[i] = fi_neq;
    }
}
