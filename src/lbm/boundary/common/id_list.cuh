#pragma once

template <auto... Ids>
struct IdList
{
    static constexpr int size = sizeof...(Ids);
};

template <class List>
struct ForEachId;

template <auto... Ids>
struct ForEachId<IdList<Ids...>>
{
    template <class F>
    __host__ __device__ static __forceinline__ void apply(F &&f)
    {
        (f.template operator()<Ids>(), ...);
    }
};

template <class List, class F>
__host__ __device__ __forceinline__ void for_each_id(F &&f)
{
    ForEachId<List>::apply(static_cast<F &&>(f));
}