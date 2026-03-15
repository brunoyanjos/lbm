// #pragma once

// #include "../../../core/types.cuh"
// #include "index_of_id.cuh"

// template <class EqList>
// struct IncomingMoments
// {
//     real_t v[EqList::size] = {};

//     template <auto Id>
//     __host__ __device__ __forceinline__ real_t &get()
//     {
//         constexpr int k = IndexOfId<Id, EqList>::value;
//         static_assert(k >= 0, "IncomingMoments::get<Id>(): Id not present in EqList.");
//         return v[k];
//     }

//     template <auto Id>
//     __host__ __device__ __forceinline__ const real_t &get() const
//     {
//         constexpr int k = IndexOfId<Id, EqList>::value;
//         static_assert(k >= 0, "IncomingMoments::get<Id>(): Id not present in EqList.");
//         return v[k];
//     }
// };