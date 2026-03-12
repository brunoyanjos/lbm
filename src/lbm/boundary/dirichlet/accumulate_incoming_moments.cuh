#pragma once

#include "../../../core/types.cuh"
#include "../common/id_list.cuh"
#include "../common/moment_basis_value.cuh"
#include "../common/incoming_moments.cuh"
#include "../common/linear_factor_value.cuh"

template <class IncomingList>
struct OpIncoming
{
    IncomingMoments<IncomingList> &incomings;
    real_t pop_i;
    const Stencil::Basis &B;

    template <auto Id>
    __device__ __forceinline__ void operator()()
    {
        incomings.template get<Id>() += pop_i * moment_basis_value(Id, B);
    }
};

template <class IncomingList>
__device__ __forceinline__ void accumulate_incoming_moments(
    IncomingMoments<IncomingList> &incomings,
    real_t pop_i,
    const Stencil::Basis &B)
{
    for_each_id<IncomingList>(OpIncoming<IncomingList>{incomings, pop_i, B});
}