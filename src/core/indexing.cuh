#pragma once

#include "geometry.h"
#include "local_domain.cuh"
#include <cstddef>
#include <limits>

constexpr size_t INVALID_INDEX = std::numeric_limits<size_t>::max();

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

__host__ __device__ __forceinline__ int wrap_global_y(int y)
{
    return wrap_y(y);
}

__host__ __device__ __forceinline__
    size_t
    idxGlobal(int x, int y)
{
    return static_cast<size_t>(x) +
           static_cast<size_t>(y) * static_cast<size_t>(NX);
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
    if ((unsigned)x >= (unsigned)NX || (unsigned)y >= (unsigned)NY)
        return INVALID_INDEX;
    return idxGlobal(x, y);
}

__host__ __device__ __forceinline__ size_t idxLocal(int x, int y_storage, int nx)
{
    return static_cast<size_t>(x) +
           static_cast<size_t>(y_storage) * static_cast<size_t>(nx);
}

__host__ __device__ __forceinline__ int localStorageYFromGlobal(int y_global,
                                                                int y_begin,
                                                                int halo)
{
    return (y_global - y_begin) + halo;
}

__host__ __device__ __forceinline__ size_t idxLocalDomainFromGlobal(int x,
                                                                    int y_global,
                                                                    const LocalDomain &domain)
{
    const int xp = wrap_x(x);
    const int yp = wrap_global_y(y_global);
    const int y_storage = localStorageYFromGlobal(yp, domain.y_begin, domain.halo);
    return idxLocal(xp, y_storage, domain.nx);
}

__device__ __forceinline__
    size_t
    idxThreadLocalInterior2D(int &x, int &y_global, const LocalDomain &domain)
{
    int y_local;
    threadXY(x, y_local);

    if ((unsigned)x >= (unsigned)domain.nx ||
        (unsigned)y_local >= (unsigned)domain.local_ny)
        return INVALID_INDEX;

    y_global = domain.y_begin + y_local;
    return idxLocal(x, y_local + domain.halo, domain.nx);
}
