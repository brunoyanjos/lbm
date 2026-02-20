#pragma once

#include "../core/types.cuh"

__host__ dim3 find_optimal_block(size_t max_lattices);
