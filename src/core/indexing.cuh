#pragma once

#include "../geometries/active_geometry.cuh"
#include <cstddef>
#include <limits>

constexpr size_t INVALID_INDEX = std::numeric_limits<size_t>::max();

__host__ __device__ __forceinline__ int wrap_x(int x)
{
    if (x < 0)
        return x + Geometry::NX;
    if (x >= Geometry::NX)
        return x - Geometry::NX;
    return x;
}

__host__ __device__ __forceinline__ int wrap_y(int y)
{
    if (y < 0)
        return y + Geometry::NY;
    if (y >= Geometry::NY)
        return y - Geometry::NY;
    return y;
}

__host__ __device__ __forceinline__
    size_t
    idxGlobal(int x, int y)
{
    return static_cast<size_t>(x) +
           static_cast<size_t>(y) * static_cast<size_t>(Geometry::NX);
}

__host__ __device__ __forceinline__
    size_t
    idxGlobalPeriodic(int x, int y)
{
    const int xp = wrap_x(x);
    const int yp = wrap_y(y);
    return idxGlobal(xp, yp);
}

__device__ __forceinline__ void threadXY(int &x, int &y)
{
    x = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
        static_cast<int>(threadIdx.x);
    y = static_cast<int>(blockIdx.y) * static_cast<int>(blockDim.y) +
        static_cast<int>(threadIdx.y);
}

__device__ __forceinline__
    size_t
    idxThreadGlobal2D(int &x, int &y)
{
    threadXY(x, y);

    if ((unsigned)x >= (unsigned)Geometry::NX || (unsigned)y >= (unsigned)Geometry::NY)
        return INVALID_INDEX;

    return idxGlobal(x, y);
}
