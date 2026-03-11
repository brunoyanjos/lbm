#include "lbm_state.cuh"
#include "../../core/cuda_utils.cuh"
#include "../../geometries/active_geometry.cuh"
#include "../../core/indexing.cuh"
#include "../../core/memory.cuh"

LBMState lbm_allocate_state()
{
    LBMState S{};
    S.N = static_cast<size_t>(Geometry::NX) * static_cast<size_t>(Geometry::NY);
    S.bytes_field = S.N * sizeof(real_t);
    S.cur = 0;

    hostMalloc_safe(S.h_rho, S.bytes_field);
    hostMalloc_safe(S.h_ux, S.bytes_field);
    hostMalloc_safe(S.h_uy, S.bytes_field);

    cudaMalloc2_safe(S.d_rho, S.bytes_field);
    cudaMalloc2_safe(S.d_ux, S.bytes_field);
    cudaMalloc2_safe(S.d_uy, S.bytes_field);

    hostMalloc_safe(S.h2.mxx, S.bytes_field);
    hostMalloc_safe(S.h2.mxy, S.bytes_field);
    hostMalloc_safe(S.h2.myy, S.bytes_field);

    cudaMalloc2_safe(S.d2.mxx, S.bytes_field);
    cudaMalloc2_safe(S.d2.mxy, S.bytes_field);
    cudaMalloc2_safe(S.d2.myy, S.bytes_field);

#if LBM_HAS_REG3_CROSS
    hostMalloc_safe(S.h3c.mxxy, S.bytes_field);
    hostMalloc_safe(S.h3c.mxyy, S.bytes_field);

    cudaMalloc2_safe(S.d3c.mxxy, S.bytes_field);
    cudaMalloc2_safe(S.d3c.mxyy, S.bytes_field);
#endif

#if LBM_HAS_REG3_AXIAL
    hostMalloc_safe(S.h3a.mxxx, S.bytes_field);
    hostMalloc_safe(S.h3a.myyy, S.bytes_field);

    cudaMalloc2_safe(S.d3a.mxxx, S.bytes_field);
    cudaMalloc2_safe(S.d3a.myyy, S.bytes_field);
#endif

    return S;
}

void lbm_free_state(LBMState &S)
{
    hostFree_safe(S.h_rho);
    hostFree_safe(S.h_ux);
    hostFree_safe(S.h_uy);

    cudaFree2_safe(S.d_rho);
    cudaFree2_safe(S.d_ux);
    cudaFree2_safe(S.d_uy);

    hostFree_safe(S.h2.mxx);
    hostFree_safe(S.h2.mxy);
    hostFree_safe(S.h2.myy);

    cudaFree2_safe(S.d2.mxx);
    cudaFree2_safe(S.d2.mxy);
    cudaFree2_safe(S.d2.myy);

#if LBM_HAS_REG3_CROSS
    hostFree_safe(S.h3c.mxxy);
    hostFree_safe(S.h3c.mxyy);

    cudaFree2_safe(S.d3c.mxxy);
    cudaFree2_safe(S.d3c.mxyy);
#endif

#if LBM_HAS_REG3_AXIAL
    hostFree_safe(S.h3a.mxxx);
    hostFree_safe(S.h3a.myyy);

    cudaFree2_safe(S.d3a.mxxx);
    cudaFree2_safe(S.d3a.myyy);
#endif

    S = {};
}