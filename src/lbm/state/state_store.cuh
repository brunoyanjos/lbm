#pragma once

#include "../../core/types.cuh"
#include "../../geometries/active_geometry.cuh"

#include "../meta/id_list.cuh"
#include "../moment/moment_id.cuh"
#include "../moment/moment_values.cuh"
#include "../state/lbm_state.cuh"

template <class MomentList>
struct NextStateStorer
{
    LBMState &S;
    int n;
    size_t idx;
    const MomentValues<MomentList> &M;

    template <auto Id>
    __device__ __forceinline__ void operator()()
    {
        if constexpr (Id == MomentId::rho)
        {
            device_field<Id>(S)[n][idx] =
                M.template get<Id>() - Geometry::RHO_0;
        }
        else
        {
            device_field<Id>(S)[n][idx] =
                M.template get<Id>();
        }
    }
};

template <class MomentList>
__device__ __forceinline__ void store_next_state(LBMState &S,
                                                 int n,
                                                 size_t idx,
                                                 const MomentValues<MomentList> &M)
{
    for_each_id<MomentList>(NextStateStorer<MomentList>{S, n, idx, M});
}