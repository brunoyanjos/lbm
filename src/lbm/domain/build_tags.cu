#include "build_tags.cuh"

#include "../../core/geometry.h"
#include "../../core/indexing.cuh"
#include "../stencil_active.cuh"
#include "../../core/cuda_utils.cuh"

#include "lbm/domain/mask_utils.cuh"

#include <cstdlib>
#include <new>

DomainTags domain_tags_allocate()
{
    DomainTags T{};
    T.N = static_cast<size_t>(NX) * static_cast<size_t>(NY);
    T.bytes_valid = T.N * sizeof(mask_t);
    T.bytes_node = T.N * sizeof(uint8_t);

    T.h_valid = static_cast<mask_t *>(std::malloc(T.bytes_valid));
    T.h_node = static_cast<uint8_t *>(std::malloc(T.bytes_node));
    if (!T.h_valid || !T.h_node)
    {
        std::free(T.h_valid);
        std::free(T.h_node);
        throw std::bad_alloc();
    }

    CUDA_CHECK(cudaMalloc(&T.d_valid, T.bytes_valid));
    CUDA_CHECK(cudaMalloc(&T.d_node, T.bytes_node));

    CUDA_CHECK(cudaMemset(T.d_valid, 0, T.bytes_valid));
    CUDA_CHECK(cudaMemset(T.d_node, 0, T.bytes_node));

    return T;
}

void domain_tags_free(DomainTags &T)
{
    std::free(T.h_valid);
    std::free(T.h_node);
    T.h_valid = nullptr;
    T.h_node = nullptr;

    if (T.d_valid)
        CUDA_CHECK(cudaFree(T.d_valid));
    if (T.d_node)
        CUDA_CHECK(cudaFree(T.d_node));
    T.d_valid = nullptr;
    T.d_node = nullptr;

    T.N = 0;
    T.bytes_valid = 0;
    T.bytes_node = 0;
}

__global__ void cavity_square_tags_kernel(mask_t *__restrict__ valid,
                                          uint8_t *__restrict__ node)
{
    int x, y;
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    const bool on_left = (x == 0);
    const bool on_right = (x == NX - 1);
    const bool on_bottom = (y == 0);
    const bool on_top = (y == NY - 1);

    const int bc_count = int(on_left) + int(on_right) + int(on_bottom) + int(on_top);

    uint8_t wid = to_u8(NodeId::FLUID);

    if (bc_count > 0)
        wid = to_u8(NodeId::DIRICHLET);

    node[idx] = wid;

    mask_t m = mask_t(0);
    m |= (mask_t(1) << 0);

#pragma unroll
    for (int i = 1; i < Stencil::Q; ++i)
    {
        const int xn = x + Stencil::cx(i);
        const int yn = y + Stencil::cy(i);

        if (xn < 0 || xn >= NX || yn < 0 || yn >= NY)
            continue;

        m |= bit(i);
    }

    valid[idx] = m;
}

void build_tags(DomainTags &T)
{
    dim3 block(16, 16, 1);
    dim3 grid((NX + block.x - 1) / block.x,
              (NY + block.y - 1) / block.y, 1);

    cavity_square_tags_kernel<<<grid, block>>>(T.d_valid, T.d_node);
    CUDA_CHECK(cudaGetLastError());

    if (T.h_valid && T.h_node)
    {
        CUDA_CHECK(cudaMemcpy(T.h_valid, T.d_valid, T.bytes_valid, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(T.h_node, T.d_node, T.bytes_node, cudaMemcpyDeviceToHost));
    }
}
