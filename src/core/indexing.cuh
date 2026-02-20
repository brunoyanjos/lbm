#pragma once

#include "geometry.h"
#include <cstddef>
#include <limits>

constexpr size_t INVALID_INDEX = std::numeric_limits<size_t>::max();

// =====================================================
// Fast periodic wrap (no modulo)
// Valid when x is in [-K, NX-1+K] with small K (e.g., K<=3).
// =====================================================

__host__ __device__ __forceinline__ int wrap_x(int x)
{
    if (x < 0)
        return x + NX;
    if (x >= NX)
        return x - NX;
    return x;
}

__host__ __device__ __forceinline__ int wrap_y(int y)
{
    if (y < 0)
        return y + NY;
    if (y >= NY)
        return y - NY;
    return y;
}

// Row-major global index (requires 0<=x<NX and 0<=y<NY)
__host__ __device__ __forceinline__
    size_t
    idxGlobal(int x, int y)
{
    return static_cast<size_t>(x) +
           static_cast<size_t>(y) * static_cast<size_t>(NX);
}

// Periodic index (safe for neighbor lookups with small offsets)
__host__ __device__ __forceinline__
    size_t
    idxGlobalPeriodic(int x, int y)
{
    const int xp = wrap_x(x);
    const int yp = wrap_y(y);
    return idxGlobal(xp, yp);
}

// Compute (x,y) from CUDA launch geometry (2D grid/block)
__device__ __forceinline__ void threadXY(int &x, int &y)
{
    x = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
        static_cast<int>(threadIdx.x);
    y = static_cast<int>(blockIdx.y) * static_cast<int>(blockDim.y) +
        static_cast<int>(threadIdx.y);
}

// Returns the global domain index for the current thread.
// Returns INVALID_INDEX if (x,y) is out of bounds.
__device__ __forceinline__
    size_t
    idxThreadGlobal2D(int &x, int &y)
{
    threadXY(x, y);
    if (x < 0 || y < 0 || x >= NX || y >= NY)
        return INVALID_INDEX;
    return idxGlobal(x, y);
}
