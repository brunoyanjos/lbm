#pragma once
#include "../../../core/types.cuh"
#include "../../../core/lbm_features.cuh"

struct Factors
{
    real_t w = r::zero;

    real_t Hx = r::zero;
    real_t Hy = r::zero;

    real_t Hxx = r::zero;
    real_t Hxy = r::zero;
    real_t Hyy = r::zero;

#if LBM_HAS_REG3_CROSS
    real_t Hxxy = r::zero;
    real_t Hxyy = r::zero;
#endif

#if LBM_HAS_REG3_AXIAL
    real_t Hxxx = r::zero;
    real_t Hyyy = r::zero;
#endif
};