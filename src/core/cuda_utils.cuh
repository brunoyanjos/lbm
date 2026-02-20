#pragma once
#include <cuda_runtime.h>
#include <iostream>
#include "types.cuh"

#define CUDA_CHECK(ans)                       \
    {                                         \
        gpuAssert((ans), __FILE__, __LINE__); \
    }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
    if (code != cudaSuccess)
    {
        std::cerr << "CUDA Error: " << cudaGetErrorString(code)
                  << " at " << file << ":" << line << std::endl;

        if (abort)
            std::exit(code);
    }
}

__device__ __forceinline__ void report_gauss_fail(int *__restrict__ err_flag,
                                                  unsigned int *__restrict__ err_count,
                                                  int idx)
{
    if (err_flag)
        atomicCAS(err_flag, 0, 1);
    if (err_count)
        atomicAdd(err_count, 1u);
}
