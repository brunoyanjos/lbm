#pragma once

#include "../../core/types.cuh"
#include "../stencil_active.cuh"
#include "../meta/for_each_id.cuh"
#include "moment_values.cuh"
#include "moment_basis_value.cuh"
#include "../layout/moment_eval.cuh"

template <class MomentList>
struct MomentAccumulator
{
    MomentValues<MomentList> &M;
    real_t pop_i;
    int icx;
    int icy;

    template <auto Id>
    __device__ __forceinline__ void operator()()
    {
        M.template get<Id>() += pop_i * moment_basis_value<Id>(icx, icy);
    }
};

__device__ __forceinline__ void evaluate_moments_from_pop(const real_t *__restrict__ pop,
                                                          MomentValues<MomentEvalList> &M)
{
#pragma unroll
    for (int k = 0; k < MomentEvalList::size; ++k)
        M.v[k] = r::zero;

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const real_t pop_i = pop[i];
        const int icx = Stencil::cx(i);
        const int icy = Stencil::cy(i);

        for_each_id<MomentEvalList>(MomentAccumulator<MomentEvalList>{M, pop_i, icx, icy});
    }
}