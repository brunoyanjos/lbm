#pragma once

#include "../../core/types.cuh"

template <class MomentList>
struct MomentValues
{
    real_t v[MomentList::size] = {};

    template <auto Id>
    __host__ __device__ __forceinline__ real_t &get()
    {
        constexpr int k = IndexOfId<Id, MomentList>::value;
        static_assert(k >= 0, "MomentValues::get<Id>(): Id not present.");
        return v[k];
    }

    template <auto Id>
    __host__ __device__ __forceinline__ const real_t &get() const
    {
        constexpr int k = IndexOfId<Id, MomentList>::value;
        static_assert(k >= 0, "MomentValues::get<Id>(): Id not present.");
        return v[k];
    }
};