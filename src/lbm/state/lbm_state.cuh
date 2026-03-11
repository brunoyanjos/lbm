#pragma once

#include "../../core/types.cuh"
#include "../../core/lbm_config.cuh"
#include "../../app/cuda_config.cuh"

struct MomentsHost2
{
    real_t *mxx = nullptr;
    real_t *mxy = nullptr;
    real_t *myy = nullptr;
};

struct MomentsDev2
{
    real_t *mxx[2] = {nullptr, nullptr};
    real_t *mxy[2] = {nullptr, nullptr};
    real_t *myy[2] = {nullptr, nullptr};
};

#if LBM_HAS_REG3_CROSS
struct MomentsHost3Cross
{
    real_t *mxxy = nullptr;
    real_t *mxyy = nullptr;
};

struct MomentsDev3Cross
{
    real_t *mxxy[2] = {nullptr, nullptr};
    real_t *mxyy[2] = {nullptr, nullptr};
};
#endif

#if LBM_HAS_REG3_AXIAL
struct MomentsHost3Axial
{
    real_t *mxxx = nullptr;
    real_t *myyy = nullptr;
};

struct MomentsDev3Axial
{
    real_t *mxxx[2] = {nullptr, nullptr};
    real_t *myyy[2] = {nullptr, nullptr};
};
#endif

struct LBMState
{
    real_t *h_rho = nullptr;
    real_t *h_ux = nullptr;
    real_t *h_uy = nullptr;

    real_t *d_rho[2] = {nullptr, nullptr};
    real_t *d_ux[2] = {nullptr, nullptr};
    real_t *d_uy[2] = {nullptr, nullptr};

    MomentsHost2 h2;
    MomentsDev2 d2;

#if LBM_HAS_REG3_CROSS
    MomentsHost3Cross h3c;
    MomentsDev3Cross d3c;
#endif

#if LBM_HAS_REG3_AXIAL
    MomentsHost3Axial h3a;
    MomentsDev3Axial d3a;
#endif

    int cur = 0;
    size_t N = 0;
    size_t bytes_field = 0;
};

[[nodiscard]] __host__ LBMState lbm_allocate_state();
__host__ void lbm_free_state(LBMState &S);