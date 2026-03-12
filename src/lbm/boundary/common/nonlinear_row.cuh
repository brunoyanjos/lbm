#pragma once

#include "../../../core/types.cuh"
#include "index_of_id.cuh"

template <class NonlinearList>
struct NonlinearRow
{
    real_t v[NonlinearList::size] = {};

    template <auto Id>
    __host__ __device__ __forceinline__ real_t &get()
    {
        constexpr int k = IndexOfId<Id, NonlinearList>::value;
        static_assert(k >= 0, "NonlinearRow::get<Id>(): Id not present in NonlinearList.");
        return v[k];
    }

    template <auto Id>
    __host__ __device__ __forceinline__ const real_t &get() const
    {
        constexpr int k = IndexOfId<Id, NonlinearList>::value;
        static_assert(k >= 0, "NonlinearRow::get<Id>(): Id not present in NonlinearList.");
        return v[k];
    }
};