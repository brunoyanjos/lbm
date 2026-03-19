#pragma once

#include <cstddef>

#include "core/types.cuh"

#include "lbm/boundary/common/system_data.cuh"
#include "lbm/boundary/dirichlet/dirichlet_accumulator.cuh"

template <std::size_t N>
__device__ __forceinline__ void build_dirichlet_system(SystemData<N> &S, real_t ux, real_t uy, const DirichletAccumulator &acc)
{
}