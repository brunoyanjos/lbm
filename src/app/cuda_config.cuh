#pragma once
#include <cuda_runtime.h>

struct CudaConfig
{
    dim3 block;
    dim3 grid;
    size_t shared_bytes;
};

[[nodiscard]] __host__ CudaConfig make_config(int shared_memory);
