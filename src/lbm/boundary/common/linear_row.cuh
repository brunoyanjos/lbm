// #pragma once

// #include "../../../core/types.cuh"
// #include "index_of_id.cuh"

// template <class VarList>
// struct LinearRow
// {
//     real_t v[VarList::size] = {};

//     template <auto Id>
//     __host__ __device__ __forceinline__ real_t &get()
//     {
//         constexpr int k = IndexOfId<Id, VarList>::value;
//         static_assert(k >= 0, "LinearRow::get<Id>(): Id not present in VarList.");
//         return v[k];
//     }

//     template <auto Id>
//     __host__ __device__ __forceinline__ const real_t &get() const
//     {
//         constexpr int k = IndexOfId<Id, VarList>::value;
//         static_assert(k >= 0, "LinearRow::get<Id>(): Id not present in VarList.");
//         return v[k];
//     }
// };