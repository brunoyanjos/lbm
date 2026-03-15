#pragma once

#include "../../core/types.cuh"

struct D2Q9
{
    static constexpr int Q = 9;

    static constexpr real_t cs2 = r_cast(1.0) / r_cast(3.0);
    static constexpr real_t as2 = r_cast(1.0) / cs2;
    static constexpr real_t as4 = as2 * as2;
    static constexpr real_t as6 = as4 * as2;

    __host__ __device__ static int cx(int i)
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

    __host__ __device__ static int cy(int i)
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

    __host__ __device__ static real_t w(int i)
    {
        switch (i)
        {
        case 0:
            return r_cast(4.0) / r_cast(9.0);
        case 1:
        case 2:
        case 3:
        case 4:
            return r_cast(1.0) / r_cast(9.0);
        default:
            return r_cast(1.0) / r_cast(36.0); // diagonais 5..8
        }
    }

    __host__ __device__ static int opp(int i)
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
};
