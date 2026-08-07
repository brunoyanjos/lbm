#pragma once
#include "lbm/moment/scale_factor.cuh"
#include "../../core/types.cuh"
#include "../../core/physics.h"

template <int RegOrder, bool HighOrder>
__device__ void moment_space_collision(NodeMomentsFor<RegOrder, HighOrder> &M)
{

    const real_t rho = M.rhoA + M.rhoB;

    const real_t inv_rho = r::one / rho;
    const real_t xA = M.rhoA * inv_rho;
    const real_t xB = M.rhoB * inv_rho;

    const real_t MASS_DENSITY_A = MASS_A * M.rhoA;
    const real_t MASS_DENSITY_B = MASS_B * M.rhoB;

    const real_t TOTAL_MASS_DENSITY = MASS_DENSITY_A + MASS_DENSITY_B;

    const real_t OMEGA_MIX = xA * OMEGA_A + xB * OMEGA_B;

    const real_t ux_col = M.ux + r::half * GRAVITY_X;
    const real_t uy_col = M.uy + r::half * GRAVITY_Y;

    const real_t mxx_col = (r::one - OMEGA_MIX) * M.mxx +
                           r::two * (r::one - OMEGA_MIX * r::half) * GRAVITY_X * M.ux;

    const real_t mxy_col = (r::one - OMEGA_MIX) * M.mxy +
                           (r::one - OMEGA_MIX * r::half) * (GRAVITY_Y * M.ux + GRAVITY_X * M.uy);

    const real_t myy_col = (r::one - OMEGA_MIX) * M.myy +
                           r::two * (r::one - OMEGA_MIX * r::half) * GRAVITY_Y * M.uy;

    M.ux = ux_col;
    M.uy = uy_col;

    M.mxx = mxx_col;
    M.mxy = mxy_col;
    M.myy = myy_col;
}
