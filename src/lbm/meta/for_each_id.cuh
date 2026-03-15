#pragma once

#include "id_list.cuh"

template <class List>
struct ForEachId;

template <auto... Ids>
struct ForEachId<IdList<Ids...>>
{
    template <class F>
    __host__ __device__ static inline void apply(F &&f)
    {
        (f.template operator()<Ids>(), ...);
    }
};

template <class List, class F>
__host__ __device__ inline void for_each_id(F &&f)
{
    ForEachId<List>::apply(static_cast<F &&>(f));
}

template <class List>
struct ForEachIdHost;

template <auto... Ids>
struct ForEachIdHost<IdList<Ids...>>
{
    template <class F>
    static inline void apply(F &&f)
    {
        (f.template operator()<Ids>(), ...);
    }
};

template <class List, class F>
inline void for_each_id_host(F &&f)
{
    ForEachIdHost<List>::apply(static_cast<F &&>(f));
}