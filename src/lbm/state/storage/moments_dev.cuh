#pragma once

#include "../../meta/index_of_id.cuh"

template <class FieldList>
struct DeviceFieldStorage
{
    real_t *v[FieldList::size][2] = {};

    template <auto Id>
    __host__ __device__ __forceinline__ real_t *(&get())[2]
    {
        constexpr int k = IndexOfId<Id, FieldList>::value;
        static_assert(k >= 0, "DeviceFieldStorage::get<Id>(): Id not present.");
        return v[k];
    }

    template <auto Id>
    __host__ __device__ __forceinline__ real_t *const (&get() const)[2]
    {
        constexpr int k = IndexOfId<Id, FieldList>::value;
        static_assert(k >= 0, "DeviceFieldStorage::get<Id>(): Id not present.");
        return v[k];
    }
};