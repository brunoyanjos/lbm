#include "lbm_state.cuh"
#include "../../core/cuda_utils.cuh"
#include "../../geometries/active_geometry.cuh"
#include "../../core/indexing.cuh"
#include "../../core/memory.cuh"
#include "../meta/for_each_id.cuh"

#include <cstdlib>
#include <cstring>
#include <new>

template <class State>
struct HostAllocator
{
    State &S;
    size_t bytes;

    template <auto Id>
    __host__ void operator()()
    {
        hostMalloc_safe(host_field<Id>(S), bytes);
    }
};

template <class State>
struct DeviceAllocator
{
    State &S;
    size_t bytes;

    template <auto Id>
    __host__ void operator()()
    {
        cudaMalloc2_safe(device_field<Id>(S), bytes);
    }
};

template <class State>
struct HostDestructor
{
    State &S;

    template <auto Id>
    __host__ void operator()()
    {
        hostFree_safe(host_field<Id>(S));
    }
};

template <class State>
struct DeviceDestructor
{
    State &S;

    template <auto Id>
    __host__ void operator()()
    {
        cudaFree2_safe(device_field<Id>(S));
    }
};

LBMState lbm_allocate_state()
{
    LBMState S{};
    S.N = static_cast<size_t>(Geometry::NX) * static_cast<size_t>(Geometry::NY);
    S.bytes_field = S.N * sizeof(real_t);
    S.cur = 0;

    try
    {
        for_each_id_host<StateFieldList>(HostAllocator<LBMState>{S, S.bytes_field});
        for_each_id_host<StateFieldList>(DeviceAllocator<LBMState>{S, S.bytes_field});
    }
    catch (...)
    {
        for_each_id_host<StateFieldList>(DeviceDestructor<LBMState>{S});
        for_each_id_host<StateFieldList>(HostDestructor<LBMState>{S});
        throw;
    }

    return S;
}

void lbm_free_state(LBMState &S)
{
    for_each_id_host<StateFieldList>(DeviceDestructor<LBMState>{S});
    for_each_id_host<StateFieldList>(HostDestructor<LBMState>{S});
    S = {};
}