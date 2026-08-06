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

    hostMalloc_safe(S.h_rhoA, S.bytes_field);
    hostMalloc_safe(S.h_uxA, S.bytes_field);
    hostMalloc_safe(S.h_uyA, S.bytes_field);
    hostMalloc_safe(S.h_mxxA, S.bytes_field);
    hostMalloc_safe(S.h_mxyA, S.bytes_field);
    hostMalloc_safe(S.h_myyA, S.bytes_field);

    hostMalloc_safe(S.h_rhoB, S.bytes_field);
    hostMalloc_safe(S.h_uxB, S.bytes_field);
    hostMalloc_safe(S.h_uyB, S.bytes_field);
    hostMalloc_safe(S.h_mxxB, S.bytes_field);
    hostMalloc_safe(S.h_mxyB, S.bytes_field);
    hostMalloc_safe(S.h_myyB, S.bytes_field);

    cudaMalloc2_safe(S.d_rhoA, S.bytes_field);
    cudaMalloc2_safe(S.d_uxA, S.bytes_field);
    cudaMalloc2_safe(S.d_uyA, S.bytes_field);
    cudaMalloc2_safe(S.d_mxxA, S.bytes_field);
    cudaMalloc2_safe(S.d_mxyA, S.bytes_field);
    cudaMalloc2_safe(S.d_myyA, S.bytes_field);

    cudaMalloc2_safe(S.d_rhoB, S.bytes_field);
    cudaMalloc2_safe(S.d_uxB, S.bytes_field);
    cudaMalloc2_safe(S.d_uyB, S.bytes_field);
    cudaMalloc2_safe(S.d_mxxB, S.bytes_field);
    cudaMalloc2_safe(S.d_mxyB, S.bytes_field);
    cudaMalloc2_safe(S.d_myyB, S.bytes_field);

    allocate_extra_state_fields(S);

    return S;
}

void lbm_free_state(LBMState &S)
{
    hostFree_safe(S.h_rhoA);
    hostFree_safe(S.h_uxA);
    hostFree_safe(S.h_uyA);
    hostFree_safe(S.h_mxxA);
    hostFree_safe(S.h_mxyA);
    hostFree_safe(S.h_myyA);

    hostFree_safe(S.h_rhoB);
    hostFree_safe(S.h_uxB);
    hostFree_safe(S.h_uyB);
    hostFree_safe(S.h_mxxB);
    hostFree_safe(S.h_mxyB);
    hostFree_safe(S.h_myyB);

    cudaFree2_safe(S.d_rhoA);
    cudaFree2_safe(S.d_uxA);
    cudaFree2_safe(S.d_uyA);
    cudaFree2_safe(S.d_mxxA);
    cudaFree2_safe(S.d_mxyA);
    cudaFree2_safe(S.d_myyA);

    cudaFree2_safe(S.d_rhoB);
    cudaFree2_safe(S.d_uxB);
    cudaFree2_safe(S.d_uyB);
    cudaFree2_safe(S.d_mxxB);
    cudaFree2_safe(S.d_mxyB);
    cudaFree2_safe(S.d_myyB);

    free_extra_state_fields(S);
}
