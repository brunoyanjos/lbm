#pragma once

#include "core/types.cuh"

namespace D2V37
{
    constexpr int Q = 37;

    constexpr real_t W0 = static_cast<real_t>(0.233150669132352502286506704066849951);
    constexpr real_t W1 = static_cast<real_t>(0.107306091542219002412464287183139936);
    constexpr real_t W2 = static_cast<real_t>(0.0576678598887948820300692153933394147);
    constexpr real_t W3 = static_cast<real_t>(0.0142082161584507502646989423441344550);
    constexpr real_t W4 = static_cast<real_t>(0.00535304900051377523273150166218849354);
    constexpr real_t W5 = static_cast<real_t>(0.00101193759267357547541090850663911471);
    constexpr real_t W6 = static_cast<real_t>(0.000245301027757717345465916643267844464);
    constexpr real_t W7 = static_cast<real_t>(0.000283414252994198217400525294194880233);

    constexpr real_t as = static_cast<real_t>(1.19697977039307435897238846385327543);
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
            return 0;
        case 11:
            return -2;
        case 12:
            return 0;
        case 13:
            return 2;
        case 14:
            return 1;
        case 15:
            return -1;
        case 16:
            return -2;
        case 17:
            return -2;
        case 18:
            return -1;
        case 19:
            return 1;
        case 20:
            return 2;
        case 21:
            return 2;
        case 22:
            return -2;
        case 23:
            return -2;
        case 24:
            return 2;
        case 25:
            return 3;
        case 26:
            return 0;
        case 27:
            return -3;
        case 28:
            return 0;
        case 29:
            return 3;
        case 30:
            return 1;
        case 31:
            return -1;
        case 32:
            return -3;
        case 33:
            return -3;
        case 34:
            return -1;
        case 35:
            return 1;
        default:
            return 3; // i == 36
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
            return 0;
        case 10:
            return 2;
        case 11:
            return 0;
        case 12:
            return -2;
        case 13:
            return 1;
        case 14:
            return 2;
        case 15:
            return 2;
        case 16:
            return 1;
        case 17:
            return -1;
        case 18:
            return -2;
        case 19:
            return -2;
        case 20:
            return -1;
        case 21:
            return 2;
        case 22:
            return 2;
        case 23:
            return -2;
        case 24:
            return -2;
        case 25:
            return 0;
        case 26:
            return 3;
        case 27:
            return 0;
        case 28:
            return -3;
        case 29:
            return 1;
        case 30:
            return 3;
        case 31:
            return 3;
        case 32:
            return 1;
        case 33:
            return -1;
        case 34:
            return -3;
        case 35:
            return -3;
        default:
            return -1; // i == 36
        }
    }

    __host__ __device__ __forceinline__ real_t w(int i)
    {
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

        case 13:
        case 14:
        case 15:
        case 16:
        case 17:
        case 18:
        case 19:
        case 20:
            return W4;

        case 21:
        case 22:
        case 23:
        case 24:
            return W5;

        case 25:
        case 26:
        case 27:
        case 28:
            return W6;

        default:
            return W7; // 29..36
        }
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
            return 17;
        case 14:
            return 18;
        case 15:
            return 19;
        case 16:
            return 20;
        case 17:
            return 13;
        case 18:
            return 14;
        case 19:
            return 15;
        case 20:
            return 16;

        case 21:
            return 23;
        case 22:
            return 24;
        case 23:
            return 21;
        case 24:
            return 22;

        case 25:
            return 27;
        case 26:
            return 28;
        case 27:
            return 25;
        case 28:
            return 26;

        case 29:
            return 33;
        case 30:
            return 34;
        case 31:
            return 35;
        case 32:
            return 36;
        case 33:
            return 29;
        case 34:
            return 30;
        case 35:
            return 31;
        default:
            return 32; // i == 36
        }
    }
}