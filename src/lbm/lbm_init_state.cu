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

__global__ void init_on_device(LBMState S)
{
    int x, y;
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    const real_t rho = RHO_0;
    const real_t ux = r::zero;
    const real_t uy = r::zero;

    real_t pop[Stencil::Q];
    equilibrium(pop, rho, ux, uy);

    const int c = S.cur;

    S.d_rho[c][idx] = rho - RHO_0;
    S.d_ux[c][idx] = ux * Stencil::as2;
    S.d_uy[c][idx] = uy * Stencil::as2;

    const real_t inv_rho = r::one / rho;

    real_t mxx = r::zero;
    real_t mxy = r::zero;
    real_t myy = r::zero;

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        mxx += pop[i] * hermite<MomentId::mxx>(i);
        mxy += pop[i] * hermite<MomentId::mxy>(i);
        myy += pop[i] * hermite<MomentId::myy>(i);
    }

    S.d_mxx[c][idx] = mxx * inv_rho * inv_scale_factor<MomentId::mxx>();
    S.d_mxy[c][idx] = mxy * inv_rho * inv_scale_factor<MomentId::mxy>();
    S.d_myy[c][idx] = myy * inv_rho * inv_scale_factor<MomentId::myy>();
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

    CUDA_CHECK(cudaMemcpy(S.h_rho, S.d_rho[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_ux, S.d_ux[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_uy, S.d_uy[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_mxx, S.d_mxx[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_mxy, S.d_mxy[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_myy, S.d_myy[c], S.bytes_field, cudaMemcpyDeviceToHost));
}
