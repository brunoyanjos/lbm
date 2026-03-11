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

    hostMalloc_safe(S.h_rho, S.bytes_field);
    hostMalloc_safe(S.h_ux, S.bytes_field);
    hostMalloc_safe(S.h_uy, S.bytes_field);
    hostMalloc_safe(S.h_mxx, S.bytes_field);
    hostMalloc_safe(S.h_mxy, S.bytes_field);
    hostMalloc_safe(S.h_myy, S.bytes_field);

    cudaMalloc_safe(S.global_rho, S.bytes_field);
    cudaMalloc_safe(S.global_ux, S.bytes_field);
    cudaMalloc_safe(S.global_uy, S.bytes_field);
    cudaMalloc_safe(S.global_mxx, S.bytes_field);
    cudaMalloc_safe(S.global_mxy, S.bytes_field);
    cudaMalloc_safe(S.global_myy, S.bytes_field);

    cudaMalloc_safe(S.d_rho, static_cast<size_t>(Geometry::NX) * 3 * sizeof(real_t));
    cudaMalloc_safe(S.d_ux, static_cast<size_t>(Geometry::NX) * 3 * sizeof(real_t));
    cudaMalloc_safe(S.d_uy, static_cast<size_t>(Geometry::NX) * 3 * sizeof(real_t));
    cudaMalloc_safe(S.d_mxx, static_cast<size_t>(Geometry::NX) * 3 * sizeof(real_t));
    cudaMalloc_safe(S.d_mxy, static_cast<size_t>(Geometry::NX) * 3 * sizeof(real_t));
    cudaMalloc_safe(S.d_myy, static_cast<size_t>(Geometry::NX) * 3 * sizeof(real_t));

    return S;
}

void lbm_free_state(LBMState &S)
{
    hostFree_safe(S.h_rho);
    hostFree_safe(S.h_ux);
    hostFree_safe(S.h_uy);
    hostFree_safe(S.h_mxx);
    hostFree_safe(S.h_mxy);
    hostFree_safe(S.h_myy);

    cudaFree_safe(S.global_rho);
    cudaFree_safe(S.global_ux);
    cudaFree_safe(S.global_uy);
    cudaFree_safe(S.global_mxx);
    cudaFree_safe(S.global_mxy);
    cudaFree_safe(S.global_myy);

    cudaFree_safe(S.d_rho);
    cudaFree_safe(S.d_ux);
    cudaFree_safe(S.d_uy);
    cudaFree_safe(S.d_mxx);
    cudaFree_safe(S.d_mxy);
    cudaFree_safe(S.d_myy);
}