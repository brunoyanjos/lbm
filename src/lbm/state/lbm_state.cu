#include "lbm_state.cuh"
#include "../../core/cuda_utils.cuh"
#include "../../core/geometry.h"
#include "../../core/indexing.cuh"
#include "../../core/memory.cuh"

template <int RegOrder, bool Rec, bool HighOrder>
void allocate_extra_state_fields(LBMStateFor<RegOrder, Rec, HighOrder> &)
{
}

template <bool HighOrder>
void allocate_extra_state_fields(LBMStateFor<3, false, HighOrder> &S)
{
    hostMalloc_safe(S.h_mxxy, S.bytes_field);
    hostMalloc_safe(S.h_mxyy, S.bytes_field);
    cudaMalloc2_safe(S.d_mxxy, S.bytes_field);
    cudaMalloc2_safe(S.d_mxyy, S.bytes_field);

    if constexpr (HighOrder)
    {
        hostMalloc_safe(S.h_mxxx, S.bytes_field);
        hostMalloc_safe(S.h_myyy, S.bytes_field);
        cudaMalloc2_safe(S.d_mxxx, S.bytes_field);
        cudaMalloc2_safe(S.d_myyy, S.bytes_field);
    }
}

template <int RegOrder, bool Rec, bool HighOrder>
void free_extra_state_fields(LBMStateFor<RegOrder, Rec, HighOrder> &)
{
}

template <bool HighOrder>
void free_extra_state_fields(LBMStateFor<3, false, HighOrder> &S)
{
    hostFree_safe(S.h_mxxy);
    hostFree_safe(S.h_mxyy);
    cudaFree2_safe(S.d_mxxy);
    cudaFree2_safe(S.d_mxyy);

    if constexpr (HighOrder)
    {
        hostFree_safe(S.h_mxxx);
        hostFree_safe(S.h_myyy);
        cudaFree2_safe(S.d_mxxx);
        cudaFree2_safe(S.d_myyy);
    }
}

LBMState lbm_allocate_state()
{
    LBMState S{};
    S.N = static_cast<size_t>(NX) * static_cast<size_t>(NY);
    S.bytes_field = S.N * sizeof(real_t);
    S.cur = 0;

    hostMalloc_safe(S.h_rho, S.bytes_field);
    hostMalloc_safe(S.h_ux, S.bytes_field);
    hostMalloc_safe(S.h_uy, S.bytes_field);
    hostMalloc_safe(S.h_mxx, S.bytes_field);
    hostMalloc_safe(S.h_mxy, S.bytes_field);
    hostMalloc_safe(S.h_myy, S.bytes_field);

    cudaMalloc2_safe(S.d_rho, S.bytes_field);
    cudaMalloc2_safe(S.d_ux, S.bytes_field);
    cudaMalloc2_safe(S.d_uy, S.bytes_field);
    cudaMalloc2_safe(S.d_mxx, S.bytes_field);
    cudaMalloc2_safe(S.d_mxy, S.bytes_field);
    cudaMalloc2_safe(S.d_myy, S.bytes_field);

    allocate_extra_state_fields(S);

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

    cudaFree2_safe(S.d_rho);
    cudaFree2_safe(S.d_ux);
    cudaFree2_safe(S.d_uy);
    cudaFree2_safe(S.d_mxx);
    cudaFree2_safe(S.d_mxy);
    cudaFree2_safe(S.d_myy);

    free_extra_state_fields(S);
}
