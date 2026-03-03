#pragma once
#include <cuda_runtime.h>

struct CudaConfig
{
    dim3 block;
    dim3 grid;
    size_t shared_bytes;
};

[[nodiscard]] __host__ CudaConfig make_config(
#if defined(LBM_STENCIL_D2V17)
    int shared_memory
#endif
);
