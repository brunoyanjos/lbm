#pragma once

#include "../../core/types.cuh"
#include "../../app/cuda_config.cuh"
#include "storage/storage.cuh"
#include "../layout/state_fields.cuh"

template <auto Id, class State>
__host__ __device__ __forceinline__ real_t *&host_field(State &S)
{
    return S.h.template get<Id>();
}

template <auto Id, class State>
__host__ __device__ __forceinline__ real_t *const &host_field(const State &S)
{
    return S.h.template get<Id>();
}

template <auto Id, class State>
__host__ __device__ __forceinline__ real_t *(&device_field(State &S))[2]
{
    return S.d.template get<Id>();
}

template <auto Id, class State>
__host__ __device__ __forceinline__ real_t *const (&device_field(const State &S))[2]
{
    return S.d.template get<Id>();
}

template <class FieldList>
struct LBMStateT
{
    HostFieldStorage<FieldList> h;
    DeviceFieldStorage<FieldList> d;

    int cur = 0;
    size_t N = 0;
    size_t bytes_field = 0;
};

using LBMState = LBMStateT<StateFieldList>;

[[nodiscard]] __host__ LBMState lbm_allocate_state();
__host__ void lbm_free_state(LBMState &S);