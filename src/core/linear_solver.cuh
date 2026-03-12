#pragma once
#include "types.cuh"
#include "math_utils.cuh"

template <int N>
struct LinearSystem
{
    real_t A[N * N]{};
    real_t b[N];
    real_t x[N];

    __device__ __forceinline__ real_t &coeff(int r, int c)
    {
        return A[r * N + c];
    }

    __device__ __forceinline__ const real_t &coeff(int r, int c) const
    {
        return A[r * N + c];
    }
};

template <int N>
__device__ __forceinline__ bool gaussianElimination(LinearSystem<N> &L)
{
    static_assert(N > 0, "N must be > 0");
    const real_t eps = (real_t)1e-20;

    // ---------- Forward elimination ----------
#pragma unroll
    for (int i = 0; i < N; ++i)
    {
        int piv = i;
        real_t maxv = dabs(L.coeff(i, i));

#pragma unroll
        for (int r = i + 1; r < N; ++r)
        {
            const real_t v = dabs(L.coeff(r, i));
            if (v > maxv)
            {
                maxv = v;
                piv = r;
            }
        }

        if (maxv <= eps)
            return false;

        if (piv != i)
        {
#pragma unroll
            for (int c = i; c < N; ++c)
            {
                const real_t tmp = L.coeff(i, c);
                L.coeff(i, c) = L.coeff(piv, c);
                L.coeff(piv, c) = tmp;
            }
            const real_t tb = L.b[i];
            L.b[i] = L.b[piv];
            L.b[piv] = tb;
        }

        const real_t inv_piv = r::one / L.coeff(i, i);

#pragma unroll
        for (int r = i + 1; r < N; ++r)
        {
            const real_t f = L.coeff(r, i) * inv_piv;
            L.coeff(r, i) = r::zero;

#pragma unroll
            for (int c = i + 1; c < N; ++c)
                L.coeff(r, c) -= f * L.coeff(i, c);

            L.b[r] -= f * L.b[i];
        }
    }

    // ---------- Back substitution ----------
#pragma unroll
    for (int i = N - 1; i >= 0; --i)
    {
        real_t s = L.b[i];

#pragma unroll
        for (int c = i + 1; c < N; ++c)
            s -= L.coeff(i, c) * L.x[c];

        const real_t diag = L.coeff(i, i);
        if (dabs(diag) <= eps)
            return false;

        L.x[i] = s / diag;
    }

    return true;
}
