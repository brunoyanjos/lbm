#pragma once
#include "types.cuh"
#include "math_utils.cuh"

template <int N>
__device__ __forceinline__ int idx_flat(int r, int c)
{
    return r * N + c;
}

template <int N>
__device__ __forceinline__ bool gaussianElimination(real_t *__restrict__ A,
                                                    real_t *__restrict__ b,
                                                    real_t *__restrict__ x)
{
    static_assert(N > 0, "N must be > 0");
    const real_t eps = (real_t)1e-20;

    // ---------- Forward elimination ----------
#pragma unroll
    for (int i = 0; i < N; ++i)
    {
        // pivot parcial na coluna i
        int piv = i;
        real_t maxv = dabs(A[idx_flat<N>(i, i)]);

#pragma unroll
        for (int r = i + 1; r < N; ++r)
        {
            const real_t v = dabs(A[idx_flat<N>(r, i)]);
            if (v > maxv)
            {
                maxv = v;
                piv = r;
            }
        }

        if (maxv <= eps)
            return false;

        // swap linhas (i <-> piv): só de i em diante
        if (piv != i)
        {
#pragma unroll
            for (int c = i; c < N; ++c)
            {
                const int ai = idx_flat<N>(i, c);
                const int ap = idx_flat<N>(piv, c);
                const real_t tmp = A[ai];
                A[ai] = A[ap];
                A[ap] = tmp;
            }
            const real_t tb = b[i];
            b[i] = b[piv];
            b[piv] = tb;
        }

        const real_t inv_piv = (real_t)1 / A[idx_flat<N>(i, i)];

        // elimina abaixo
#pragma unroll
        for (int r = i + 1; r < N; ++r)
        {
            const real_t f = A[idx_flat<N>(r, i)] * inv_piv;
            A[idx_flat<N>(r, i)] = (real_t)0;

#pragma unroll
            for (int c = i + 1; c < N; ++c)
                A[idx_flat<N>(r, c)] -= f * A[idx_flat<N>(i, c)];

            b[r] -= f * b[i];
        }
    }

    // ---------- Back substitution ----------
#pragma unroll
    for (int i = N - 1; i >= 0; --i)
    {
        real_t s = b[i];

#pragma unroll
        for (int c = i + 1; c < N; ++c)
            s -= A[idx_flat<N>(i, c)] * x[c];

        const real_t diag = A[idx_flat<N>(i, i)];
        if (dabs(diag) <= eps)
            return false;

        x[i] = s / diag;
    }

    return true;
}
