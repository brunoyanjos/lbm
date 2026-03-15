#pragma once

#include "../../meta/index_of_id.cuh"

template <class FieldList>
struct HostFieldStorage
{
    real_t *v[FieldList::size] = {};

    template <auto Id>
    __host__ __device__ __forceinline__ real_t *&get()
    {
        constexpr int k = IndexOfId<Id, FieldList>::value;
        static_assert(k >= 0, "HostFieldStorage::get<Id>(): Id not present.");
        return v[k];
    }

    template <auto Id>
    __host__ __device__ __forceinline__ real_t *const &get() const
    {
        constexpr int k = IndexOfId<Id, FieldList>::value;
        static_assert(k >= 0, "HostFieldStorage::get<Id>(): Id not present.");
        return v[k];
    }
};