#pragma once
#include <cstddef>

#include "core/types.cuh"

#ifdef LBM_GRID
constexpr int NX = LBM_GRID;
constexpr int NY = LBM_GRID;
#else
constexpr int NX = 512;
constexpr int NY = 512;
#endif

constexpr real_t xc = NX / 2;
constexpr real_t yc = NY / 2;
