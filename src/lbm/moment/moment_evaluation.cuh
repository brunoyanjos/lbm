#pragma once

#include "core/types.cuh"
#include "lbm/stencil_active.cuh"
#include "lbm/hermite/hermite.cuh"
#include "lbm/moment/moment_id.cuh"
#include "lbm/moment/node_moments.cuh"

template <int Order, bool HighOrder>
__device__ __forceinline__ void evaluate_moments_from_pop(
    const real_t *__restrict__ pop, const real_t *__restrict__ popA,
    const real_t *__restrict__ popB,
    NodeMomentsFor<Order, HighOrder> &M)
{
    M.rhoA = r::zero;
    M.rhoB = r::zero;

    M.ux = r::zero;
    M.uy = r::zero;
    M.mxx = r::zero;
    M.mxy = r::zero;
    M.myy = r::zero;

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const real_t fiA = popA[i] - Stencil::w(i);
        const real_t fiB = popB[i] - Stencil::w(i);
        const real_t fi_neq = pop[i] - Stencil::w(i);

        M.rhoA += fiA;
        M.rhoB += fiB;

        M.ux += (MASS_A * fiA + MASS_B * fiB) * hermite<MomentId::ux>(i);
        M.uy += (MASS_A * fiA + MASS_B * fiB) * hermite<MomentId::uy>(i);

        M.mxx += (MASS_A * fiA + MASS_B * fiB + fi_neq) * hermite<MomentId::mxx>(i);
        M.mxy += (MASS_A * fiA + MASS_B * fiB + fi_neq) * hermite<MomentId::mxy>(i);
        M.myy += (MASS_A * fiA + MASS_B * fiB + fi_neq) * hermite<MomentId::myy>(i);
    }

    M.rhoA += RHOA_0;
    M.rhoB += RHOB_0;

    const real_t MASS_DENSITY_A = MASS_A * M.rhoA;
    const real_t MASS_DENSITY_B = MASS_B * M.rhoB;

    const real_t TOTAL_MASS_DENSITY = MASS_DENSITY_A + MASS_DENSITY_B;

    const real_t inv_mass = r::one / (MASS_DENSITY_A + MASS_DENSITY_B);
    const real_t inv_rho = r::one / (M.rhoA + M.rhoB);

    const real_t rho = M.rhoA + M.rhoB;

    const real_t THETA_A = (INV_MASS_A - r::one) * Stencil::cs2;
    const real_t THETA_B = (INV_MASS_B - r::one) * Stencil::cs2;

    M.ux *= inv_mass;
    M.uy *= inv_mass;

    M.mxx -= TOTAL_MASS_DENSITY * M.ux * M.ux + (rho - TOTAL_MASS_DENSITY) * Stencil::cs2;
    M.mxy -= TOTAL_MASS_DENSITY * M.ux * M.uy;
    M.myy -= TOTAL_MASS_DENSITY * M.uy * M.uy + (rho - TOTAL_MASS_DENSITY) * Stencil::cs2;

    M.mxx *= inv_mass;
    M.mxy *= inv_mass;
    M.myy *= inv_mass;

    M.ux += GRAVITY_X * r::half;
    M.uy += GRAVITY_Y * r::half;

    M.mxx += TOTAL_MASS_DENSITY * inv_rho * GRAVITY_X * M.ux;
    M.mxy += TOTAL_MASS_DENSITY * inv_rho * r::half * (GRAVITY_Y * M.ux + GRAVITY_X * M.uy);
    M.myy += TOTAL_MASS_DENSITY * inv_rho * GRAVITY_Y * M.uy;
}