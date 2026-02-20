#pragma once
#include <cuda_runtime.h>
#include "../core/cuda_utils.cuh"

namespace app
{

    struct BenchmarkResult
    {
        double gpu_seconds = 0.0;  // tempo medido por cudaEvent (kernel+device work enfileirado)
        double wall_seconds = 0.0; // tempo real (inclui host/IO se habilitado)
        double mlups_gpu = 0.0;
        double mlups_wall = 0.0;
        int measured_steps = 0;
    };

    class GpuTimer
    {
    public:
        void start(cudaStream_t stream = 0)
        {
            CUDA_CHECK(cudaEventCreate(&ev0));
            CUDA_CHECK(cudaEventCreate(&ev1));
            CUDA_CHECK(cudaEventRecord(ev0, stream));
            stream_ = stream;
        }

        double stop_seconds()
        {
            CUDA_CHECK(cudaEventRecord(ev1, stream_));
            CUDA_CHECK(cudaEventSynchronize(ev1));
            float ms = 0.0f;
            CUDA_CHECK(cudaEventElapsedTime(&ms, ev0, ev1));
            CUDA_CHECK(cudaEventDestroy(ev0));
            CUDA_CHECK(cudaEventDestroy(ev1));
            ev0 = ev1 = nullptr;
            return double(ms) * 1e-3;
        }

    private:
        cudaEvent_t ev0 = nullptr;
        cudaEvent_t ev1 = nullptr;
        cudaStream_t stream_ = 0;
    };

}
