#include "lbm_init_state.cuh"
#include <cuda_runtime.h>

#include "../geometries/active_geometry.cuh"
#include "../core/indexing.cuh"
#include "../core/cuda_utils.cuh"

#include "lbm_equilibrium.cuh"
#include "stencil_active.cuh"

__global__ void init_on_device(LBMState S)
{
    int x, y;
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    const real_t rho = Geometry::RHO_0;
    const real_t ux = real_t(0);
    const real_t uy = real_t(0);

    real_t pop[Stencil::Q];
    equilibrium(pop, rho, ux, uy);

    S.global_rho[idx] = rho - Geometry::RHO_0;
    S.global_ux[idx] = ux * Stencil::as2;
    S.global_uy[idx] = uy * Stencil::as2;

    const real_t inv_rho = real_t(1) / rho;

    real_t mxx = real_t(0);
    real_t mxy = real_t(0);
    real_t myy = real_t(0);

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const real_t cx = static_cast<real_t>(Stencil::cx(i));
        const real_t cy = static_cast<real_t>(Stencil::cy(i));

        const real_t Hxx = cx * cx - Stencil::cs2;
        const real_t Hxy = cx * cy;
        const real_t Hyy = cy * cy - Stencil::cs2;

        mxx += pop[i] * Hxx;
        mxy += pop[i] * Hxy;
        myy += pop[i] * Hyy;
    }

    S.global_mxx[idx] = mxx * inv_rho * (Stencil::as4 * real_t(0.5));
    S.global_mxy[idx] = mxy * inv_rho * Stencil::as4;
    S.global_myy[idx] = myy * inv_rho * (Stencil::as4 * real_t(0.5));
}

void init_state(LBMState &S, const CudaConfig &cfg)
{
    init_on_device<<<cfg.grid, cfg.block>>>(S);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

void upload_state_to_host(LBMState &S)
{
    CUDA_CHECK(cudaMemcpy(S.h_rho, S.global_rho, S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_ux, S.global_ux, S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_uy, S.global_uy, S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_mxx, S.global_mxx, S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_mxy, S.global_mxy, S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h_myy, S.global_myy, S.bytes_field, cudaMemcpyDeviceToHost));
}