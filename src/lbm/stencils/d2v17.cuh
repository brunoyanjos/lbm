#pragma once

#include "../../core/types.cuh"

namespace D2V17
{
    constexpr int Q = 17;

    constexpr real_t W0 = static_cast<real_t>(0.40200514690911262594166439245907543027);
    constexpr real_t W1 = static_cast<real_t>(0.11615486649778154387403545662591119451);
    constexpr real_t W2 = static_cast<real_t>(0.03300635362298691394975382591472168169);
    constexpr real_t W3 = static_cast<real_t>(0.00007907860216591813123713324054816357);
    constexpr real_t W4 = static_cast<real_t>(0.00025841454978746755955748610405009882);

    constexpr real_t as = static_cast<real_t>(1.64343060879795414682265416738949126364);
    constexpr real_t as2 = as * as;
    constexpr real_t as4 = as2 * as2;
    constexpr real_t as6 = as4 * as2;

    constexpr real_t cs2 = static_cast<real_t>(1.0) / as2;

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
        case 8:
            return 1;
        case 9:
            return 2;
        case 10:
            return -2;
        case 11:
            return -2;
        case 12:
            return 2;
        case 13:
            return 3;
        case 14:
            return 0;
        case 15:
            return -3;
        default:
            return 0; // i==16
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
        case 8:
            return -1;
        case 9:
            return 2;
        case 10:
            return 2;
        case 11:
            return -2;
        case 12:
            return -2;
        case 13:
            return 0;
        case 14:
            return 3;
        case 15:
            return 0;
        default:
            return -3; // i==16
        }
    }

    __host__ __device__ __forceinline__ real_t w(int i)
    {
        // agrupar pesos iguais ajuda o compilador
        switch (i)
        {
        case 0:
            return W0;
        case 1:
        case 2:
        case 3:
        case 4:
            return W1;
        case 5:
        case 6:
        case 7:
        case 8:
            return W2;
        case 9:
        case 10:
        case 11:
        case 12:
            return W3;
        default:
            return W4; // 13..16
        }
    }

    __host__ __device__ __forceinline__ void basis2(int i, real_t &cx, real_t &cy, real_t &Hxx, real_t &Hxy, real_t &Hyy)
    {
        cx = static_cast<real_t>(D2V17::cx(i));
        cy = static_cast<real_t>(D2V17::cy(i));
        Hxx = cx * cx - D2V17::cs2;
        Hxy = cx * cy;
        Hyy = cy * cy - D2V17::cs2;
    }

    __host__ __device__ __forceinline__ void basis2polar(int i, real_t c, real_t s,
                                                         real_t &cx, real_t &cy, real_t &Hxx, real_t &Hxy, real_t &Hyy)
    {
        cx = r_cast(D2V17::cx(i)) * c + r_cast(D2V17::cy(i)) * s;
        cy = r_cast(D2V17::cy(i)) * c - r_cast(D2V17::cx(i)) * s;
        Hxx = cx * cx - D2V17::cs2;
        Hxy = cx * cy;
        Hyy = cy * cy - D2V17::cs2;
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
        case 8:
            return 6;

        case 9:
            return 11;
        case 10:
            return 12;
        case 11:
            return 9;
        case 12:
            return 10;

        case 13:
            return 15;
        case 14:
            return 16;
        case 15:
            return 13;
        default:
            return 14; // i==16
        }
    }
} // namespace D2V17
