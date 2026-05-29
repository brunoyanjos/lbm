#include "lbm_state.cuh"
#include "core/cuda_utils.cuh"
#include "core/geometry.h"
#include "core/indexing.cuh"
#include "core/memory.cuh"

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

LBMState lbm_allocate_state(const LocalDomain &domain)
{
    LBMState S{};
    S.domain = domain;
    S.N = domain.nx * domain.storage_ny;
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

[[nodiscard]] __host__ std::vector<LBMState> allocate_partition_states(const app::RunContext &ctx)
{
    std::vector<LBMState> states;
    states.reserve(ctx.partitions.size());

    std::cout << "[MULTI_GPU] allocating " << ctx.partitions.size()
              << " local state partitions\n";

    for (const auto &partition : ctx.partitions)
    {
        CUDA_CHECK(cudaSetDevice(partition.device_id));
        const LocalDomain domain = make_local_domain(partition.y_begin,
                                                     partition.y_end,
                                                     partition.halo);
        states.push_back(lbm_allocate_state(domain));

        std::cout << "[MULTI_GPU] device=" << partition.device_id
                  << " y=[" << domain.y_begin << "," << domain.y_end << ")"
                  << " local_ny=" << domain.local_ny
                  << " halo=" << domain.halo
                  << " storage_ny=" << domain.storage_ny
                  << " bytes/field=" << states.back().bytes_field << "\n";
    }

    return states;
}

__host__ void free_partition_states(std::vector<LBMState> &states, const app::RunContext &ctx)
{
    for (size_t i = 0; i < states.size(); ++i)
    {
        CUDA_CHECK(cudaSetDevice(ctx.partitions[i].device_id));
        lbm_free_state(states[i]);
    }

    states.clear();
}