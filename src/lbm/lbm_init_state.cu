#include "lbm_init_state.cuh"
#include <cuda_runtime.h>

#include "core/geometry.h"
#include "core/physics.h"
#include "core/indexing.cuh"
#include "core/cuda_utils.cuh"

#include "lbm/hermite/hermite.cuh"
#include "lbm/moment/scale_factor.cuh"

#include "lbm_equilibrium.cuh"
#include "stencil_active.cuh"

template <int RegOrder, bool Rec, bool HighOrder>
__device__ __forceinline__ void init_extra_state_fields(LBMStateFor<RegOrder, Rec, HighOrder> &, int, size_t)
{
}

template <bool HighOrder>
__device__ __forceinline__ void init_extra_state_fields(LBMStateFor<3, false, HighOrder> &S, int c, size_t idx)
{
    S.d_mxxy[c][idx] = r::zero;
    S.d_mxyy[c][idx] = r::zero;

    if constexpr (HighOrder)
    {
        S.d_mxxx[c][idx] = r::zero;
        S.d_myyy[c][idx] = r::zero;
    }
}

template <int RegOrder, bool Rec, bool HighOrder>
__host__ void upload_extra_state_fields(LBMStateFor<RegOrder, Rec, HighOrder> &, int)
{
}

template <bool HighOrder>
__host__ void upload_extra_state_fields(LBMStateFor<3, false, HighOrder> &S, int c)
{
    CUDA_CHECK(cudaMemcpy(S.h_mxxy, S.d_mxxy[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_mxyy, S.d_mxyy[c], S.bytes_field, cudaMemcpyDeviceToHost));

    if constexpr (HighOrder)
    {
        CUDA_CHECK(cudaMemcpy(S.h_mxxx, S.d_mxxx[c], S.bytes_field, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(S.h_myyy, S.d_myyy[c], S.bytes_field, cudaMemcpyDeviceToHost));
    }
}

__global__ void init_on_device(LBMState S)
{
    int x, y;
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    real_t rhoA = RHOA_0;
    real_t rhoB = RHOB_0;
    const real_t ux = r::zero;
    const real_t uy = r::zero;

    // const real_t x_diff = x - NX / 2;
    // const real_t y_diff = y - NY / 2;

    // if (x_diff * x_diff + y_diff * y_diff < RADIUS * RADIUS)
    if (x < NY / 2)
    {
        rhoA = r::zero;
        // rhoB = r::zero;
    }
    else
    {
        // rhoA = r::zero;
        rhoB = r::zero;
    }

    real_t popA[Stencil::Q];
    equilibrium(popA, rhoA, ux, uy);

    real_t popB[Stencil::Q];
    equilibrium(popB, rhoB, ux, uy);

    const int c = S.cur;

    S.d_rhoA[c][idx] = rhoA;
    S.d_uxA[c][idx] = ux * Stencil::as2;
    S.d_uyA[c][idx] = uy * Stencil::as2;

    S.d_rhoB[c][idx] = rhoB;
    S.d_uxB[c][idx] = ux * Stencil::as2;
    S.d_uyB[c][idx] = uy * Stencil::as2;

    const real_t inv_rho = r::one / (rhoA + rhoB);

    real_t mxx = r::zero;
    real_t mxy = r::zero;
    real_t myy = r::zero;

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        mxx += (popA[i] + popB[i]) * hermite<MomentId::mxx>(i);
        mxy += (popA[i] + popB[i]) * hermite<MomentId::mxy>(i);
        myy += (popA[i] + popB[i]) * hermite<MomentId::myy>(i);
    }

    S.d_mxxA[c][idx] = mxx * inv_rho * scale_factor<MomentId::mxx>();
    S.d_mxyA[c][idx] = mxy * inv_rho * scale_factor<MomentId::mxy>();
    S.d_myyA[c][idx] = myy * inv_rho * scale_factor<MomentId::myy>();

    S.d_mxxB[c][idx] = mxx * inv_rho * scale_factor<MomentId::mxx>();
    S.d_mxyB[c][idx] = mxy * inv_rho * scale_factor<MomentId::mxy>();
    S.d_myyB[c][idx] = myy * inv_rho * scale_factor<MomentId::myy>();

    init_extra_state_fields(S, c, idx);
}

void init_state(LBMState &S, const CudaConfig &cfg)
{
    S.cur = 0;

    init_on_device<<<cfg.grid, cfg.block>>>(S);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

void upload_state_to_host(LBMState &S)
{
    const int c = S.cur;

    CUDA_CHECK(cudaMemcpy(S.h_rhoA, S.d_rhoA[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_uxA, S.d_uxA[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_uyA, S.d_uyA[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_mxxA, S.d_mxxA[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_mxyA, S.d_mxyA[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_myyA, S.d_myyA[c], S.bytes_field, cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaMemcpy(S.h_rhoB, S.d_rhoB[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_uxB, S.d_uxB[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_uyB, S.d_uyB[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_mxxB, S.d_mxxB[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_mxyB, S.d_mxyB[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_myyB, S.d_myyB[c], S.bytes_field, cudaMemcpyDeviceToHost));

    upload_extra_state_fields(S, c);
}
