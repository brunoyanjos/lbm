#include "lbm_init_state.cuh"
#include <cuda_runtime.h>

#include "../geometries/active_geometry.cuh"
#include "../core/indexing.cuh"
#include "../core/cuda_utils.cuh"

#include "lbm_equilibrium.cuh"
#include "stencil_active.cuh"
#include "../core/lbm_config.cuh"

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

    const int c = S.cur;

    S.d_rho[c][idx] = rho - Geometry::RHO_0;
    S.d_ux[c][idx] = ux * Stencil::as2;
    S.d_uy[c][idx] = uy * Stencil::as2;

    const real_t inv_rho = real_t(1) / rho;

    real_t mxx = real_t(0);
    real_t mxy = real_t(0);
    real_t myy = real_t(0);

#if LBM_HAS_REG3_CROSS
    real_t mxxy = real_t(0);
    real_t mxyy = real_t(0);
#endif

#if LBM_HAS_REG3_AXIAL
    real_t mxxx = real_t(0);
    real_t myyy = real_t(0);
#endif

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const auto B = Stencil::basis(i);

        mxx += pop[i] * B.Hxx;
        mxy += pop[i] * B.Hxy;
        myy += pop[i] * B.Hyy;

#if LBM_HAS_REG3_CROSS
        mxxy += pop[i] * B.Hxxy;
        mxyy += pop[i] * B.Hxyy;
#endif

#if LBM_HAS_REG3_AXIAL
        mxxx += pop[i] * B.Hxxx;
        myyy += pop[i] * B.Hyyy;
#endif
    }

    S.d2.mxx[c][idx] = mxx * inv_rho * (Stencil::as4 * r::half);
    S.d2.mxy[c][idx] = mxy * inv_rho * Stencil::as4;
    S.d2.myy[c][idx] = myy * inv_rho * (Stencil::as4 * r::half);

#if LBM_HAS_REG3_CROSS
    S.d3c.mxxy[c][idx] = mxxy * inv_rho * (Stencil::as6 * r::half);
    S.d3c.mxyy[c][idx] = mxyy * inv_rho * (Stencil::as6 * r::half);
#endif

#if LBM_HAS_REG3_AXIAL
    S.d3a.mxxx[c][idx] = mxxx * inv_rho * (Stencil::as6 * r::sixth);
    S.d3a.myyy[c][idx] = myyy * inv_rho * (Stencil::as6 * r::sixth);
#endif
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

    CUDA_CHECK(cudaMemcpy(S.h2.mxx, S.d2.mxx[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h2.mxy, S.d2.mxy[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h2.myy, S.d2.myy[c], S.bytes_field, cudaMemcpyDeviceToHost));

#if LBM_HAS_REG3_CROSS
    CUDA_CHECK(cudaMemcpy(S.h3c.mxxy, S.d3c.mxxy[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h3c.mxyy, S.d3c.mxyy[c], S.bytes_field, cudaMemcpyDeviceToHost));
#endif

#if LBM_HAS_REG3_AXIAL
    CUDA_CHECK(cudaMemcpy(S.h3a.mxxx, S.d3a.mxxx[c], S.bytes_field, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(S.h3a.myyy, S.d3a.myyy[c], S.bytes_field, cudaMemcpyDeviceToHost));
#endif
}