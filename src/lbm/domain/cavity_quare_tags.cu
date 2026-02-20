#include "cavity_square_tags.cuh"

#include "../../core/geometry.h"
#include "../../core/indexing.cuh"
#include "../stencil_active.cuh"
#include "../../core/cuda_utils.cuh"

#include <cstdlib>
#include <new>

DomainTags domain_tags_allocate(bool host_buffers)
{
    DomainTags T{};
    T.N = static_cast<size_t>(NX) * static_cast<size_t>(NY);
    T.bytes_valid = T.N * sizeof(uint32_t);
    T.bytes_wall = T.N * sizeof(uint8_t);

    if (host_buffers)
    {
        T.h_valid = static_cast<uint32_t *>(std::malloc(T.bytes_valid));
        T.h_wall = static_cast<uint8_t *>(std::malloc(T.bytes_wall));
        if (!T.h_valid || !T.h_wall)
        {
            std::free(T.h_valid);
            std::free(T.h_wall);
            throw std::bad_alloc();
        }
    }

    CUDA_CHECK(cudaMalloc(&T.d_valid, T.bytes_valid));
    CUDA_CHECK(cudaMalloc(&T.d_wall, T.bytes_wall));

    CUDA_CHECK(cudaMemset(T.d_valid, 0, T.bytes_valid));
    CUDA_CHECK(cudaMemset(T.d_wall, 0, T.bytes_wall));

    return T;
}

void domain_tags_free(DomainTags &T)
{
    std::free(T.h_valid);
    std::free(T.h_wall);
    T.h_valid = nullptr;
    T.h_wall = nullptr;

    if (T.d_valid)
        CUDA_CHECK(cudaFree(T.d_valid));
    if (T.d_wall)
        CUDA_CHECK(cudaFree(T.d_wall));
    T.d_valid = nullptr;
    T.d_wall = nullptr;

    T.N = 0;
    T.bytes_valid = 0;
    T.bytes_wall = 0;
}

__global__ void cavity_square_tags_kernel(uint32_t *__restrict__ valid,
                                          uint8_t *__restrict__ wall)
{
    int x, y;
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    // wall-id (para Dirichlet: TOP tem tampa móvel; demais u=0)
    const bool on_left = (x == 0);
    const bool on_right = (x == NX - 1);
    const bool on_bottom = (y == 0);
    const bool on_top = (y == NY - 1);

    const int bc_count = int(on_left) + int(on_right) + int(on_bottom) + int(on_top);

    uint8_t wid = static_cast<uint8_t>(WallId::NONE);

    if (bc_count >= 2)
        wid = static_cast<uint8_t>(WallId::CORNER);
    else if (on_left)
        wid = static_cast<uint8_t>(WallId::LEFT);
    else if (on_right)
        wid = static_cast<uint8_t>(WallId::RIGHT);
    else if (on_bottom)
        wid = static_cast<uint8_t>(WallId::BOTTOM);
    else if (on_top)
        wid = static_cast<uint8_t>(WallId::TOP);

    wall[idx] = wid;

    // valid directions: vizinho dentro do domínio
    uint32_t m = 0u;
    m |= (1u << 0); // rest sempre válido

#pragma unroll
    for (int i = 1; i < Stencil::Q; ++i)
    {
        const int xn = x + Stencil::cx(i);
        const int yn = y + Stencil::cy(i);

        if (xn < 0 || xn >= NX || yn < 0 || yn >= NY)
            continue;

        m |= (1u << i);
    }

    valid[idx] = m;
}

void build_cavity_square_tags(DomainTags &T)
{
    dim3 block(16, 16, 1);
    dim3 grid((NX + block.x - 1) / block.x,
              (NY + block.y - 1) / block.y, 1);

    cavity_square_tags_kernel<<<grid, block>>>(T.d_valid, T.d_wall);
    CUDA_CHECK(cudaGetLastError());

    if (T.h_valid && T.h_wall)
    {
        CUDA_CHECK(cudaMemcpy(T.h_valid, T.d_valid, T.bytes_valid, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(T.h_wall, T.d_wall, T.bytes_wall, cudaMemcpyDeviceToHost));
    }
}
