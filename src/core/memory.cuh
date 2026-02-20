#pragma once
#include <cuda_runtime.h>
#include <cstddef>

template <typename T>
__host__ inline void cudaMalloc2_safe(T *(&p)[2], size_t bytes)
{
    p[0] = nullptr;
    p[1] = nullptr;

    cudaError_t e0 = cudaMalloc((void **)&p[0], bytes);
    if (e0 != cudaSuccess)
    {
        p[0] = nullptr;
        p[1] = nullptr;
        CUDA_CHECK(e0);
    }

    cudaError_t e1 = cudaMalloc((void **)&p[1], bytes);
    if (e1 != cudaSuccess)
    {
        cudaFree(p[0]);
        p[0] = nullptr;
        p[1] = nullptr;
        CUDA_CHECK(e1);
    }

    cudaError_t m0 = cudaMemset(p[0], 0, bytes);
    if (m0 != cudaSuccess)
    {
        cudaFree(p[1]);
        cudaFree(p[0]);
        p[0] = nullptr;
        p[1] = nullptr;
        CUDA_CHECK(m0);
    }

    cudaError_t m1 = cudaMemset(p[1], 0, bytes);
    if (m1 != cudaSuccess)
    {
        cudaFree(p[1]);
        cudaFree(p[0]);
        p[0] = nullptr;
        p[1] = nullptr;
        CUDA_CHECK(m1);
    }
}

template <typename T>
__host__ inline void cudaFree2_safe(T *(&p)[2])
{
    if (p[0])
        CUDA_CHECK(cudaFree(p[0]));
    if (p[1])
        CUDA_CHECK(cudaFree(p[1]));
    p[0] = nullptr;
    p[1] = nullptr;
}

template <typename T>
__host__ inline void hostMalloc_safe(T *&p, size_t bytes)
{
    p = static_cast<T *>(std::malloc(bytes));
    if (!p)
        throw std::bad_alloc();
}

template <typename T>
__host__ inline void hostFree_safe(T *&p)
{
    std::free(p);
    p = nullptr;
}
