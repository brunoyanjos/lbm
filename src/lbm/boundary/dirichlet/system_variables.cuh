#pragma once
#include "../../../core/types.cuh"

struct DirichletVariable2D
{
    // incoming raw moments
    real_t rho_I = r::zero;

    real_t mxx_I = r::zero;
    real_t mxy_I = r::zero;
    real_t myy_I = r::zero;

    // outgoing sums from rho equation
    real_t rho_rho = r::zero;
    real_t rho_ux = r::zero;
    real_t rho_uy = r::zero;
    real_t rho_mxx = r::zero;
    real_t rho_mxy = r::zero;
    real_t rho_myy = r::zero;
    real_t rho_uxux = r::zero;
    real_t rho_uxuy = r::zero;
    real_t rho_uyuy = r::zero;

    // mxx_I equation
    real_t mxx_rho = r::zero;
    real_t mxx_ux = r::zero;
    real_t mxx_uy = r::zero;
    real_t mxx_mxx = r::zero;
    real_t mxx_mxy = r::zero;
    real_t mxx_myy = r::zero;

    // mxy_I equation
    real_t mxy_rho = r::zero;
    real_t mxy_ux = r::zero;
    real_t mxy_uy = r::zero;
    real_t mxy_mxx = r::zero;
    real_t mxy_mxy = r::zero;
    real_t mxy_myy = r::zero;

    // myy_I equation
    real_t myy_rho = r::zero;
    real_t myy_ux = r::zero;
    real_t myy_uy = r::zero;
    real_t myy_mxx = r::zero;
    real_t myy_mxy = r::zero;
    real_t myy_myy = r::zero;

    __device__ __forceinline__ void normalize_moments()
    {
        const real_t inv_rho = r::one / rho_I;

        mxx_I *= inv_rho;
        mxy_I *= inv_rho;
        myy_I *= inv_rho;
    }
};
