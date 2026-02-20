#pragma once

#include "../../core/types.cuh"

namespace D2Q9
{
    constexpr int Q = 9;

    // Para D2Q9 clássico (cs^2 = 1/3)
    constexpr real_t cs2 = static_cast<real_t>(1.0) / static_cast<real_t>(3.0);
    constexpr real_t as2 = static_cast<real_t>(1.0) / cs2; // = 3
    constexpr real_t as4 = as2 * as2;

    __host__ __device__ __forceinline__ int cx(int i)
    {
        switch (i)
        {
        case 0:
            return 0;
        case 1:
            return 1;
        case 2:
            return 0;
        case 3:
            return -1;
        case 4:
            return 0;
        case 5:
            return 1;
        case 6:
            return -1;
        case 7:
            return -1;
        default:
            return 1; // i==8
        }
    }

    __host__ __device__ __forceinline__ int cy(int i)
    {
        switch (i)
        {
        case 0:
            return 0;
        case 1:
            return 0;
        case 2:
            return 1;
        case 3:
            return 0;
        case 4:
            return -1;
        case 5:
            return 1;
        case 6:
            return 1;
        case 7:
            return -1;
        default:
            return -1; // i==8
        }
    }

    __host__ __device__ __forceinline__ real_t w(int i)
    {
        switch (i)
        {
        case 0:
            return static_cast<real_t>(4.0) / static_cast<real_t>(9.0);
        case 1:
        case 2:
        case 3:
        case 4:
            return static_cast<real_t>(1.0) / static_cast<real_t>(9.0);
        default:
            return static_cast<real_t>(1.0) / static_cast<real_t>(36.0); // diagonais 5..8
        }
    }

    __host__ __device__ __forceinline__ void basis2(int i, real_t &cx, real_t &cy, real_t &Hxx, real_t &Hxy, real_t &Hyy)
    {
        cx = static_cast<real_t>(D2Q9::cx(i));
        cy = static_cast<real_t>(D2Q9::cy(i));
        Hxx = cx * cx - D2Q9::cs2;
        Hxy = cx * cy;
        Hyy = cy * cy - D2Q9::cs2;
    }

    __host__ __device__ __forceinline__ int opp(int i)
    {
        switch (i)
        {
        case 0:
            return 0;
        case 1:
            return 3;
        case 2:
            return 4;
        case 3:
            return 1;
        case 4:
            return 2;
        case 5:
            return 7;
        case 6:
            return 8;
        case 7:
            return 5;
        default:
            return 6; // i==8
        }
    }
}
