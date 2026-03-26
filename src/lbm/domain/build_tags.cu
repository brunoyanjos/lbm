#include "build_tags.cuh"

#include "core/types.cuh"
#include "core/geometry.h"
#include "core/physics.h"
#include "core/indexing.cuh"
#include "core/math_utils.cuh"
#include "../stencil_active.cuh"
#include "core/cuda_utils.cuh"

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

__global__ void annul_tags_kernel(uint8_t *__restrict__ nodes)
{
    int x, y;
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    const real_t x_center = r_cast(NX - 1) * r::half;
    const real_t y_center = r_cast(NY - 1) * r::half;

    const real_t dx = x - x_center;
    const real_t dy = y - y_center;

    const real_t radius = r_sqrt(dx * dx + dy * dy);

    nodes[idx] = to_u8(NodeId::FLUID);

    if (radius < R_IN || radius > R_OUT)
        nodes[idx] = to_u8(NodeId::SOLID);
}

__global__ void annul_tags_boundary_kernel(uint8_t *__restrict__ nodes)
{
    int x, y;
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    const uint8_t node_id = nodes[idx];

    if (node_id == to_u8(NodeId::FLUID))
        return;

#pragma unroll
    for (int i = 1; i < Stencil::Q; ++i)
    {
        uint8_t adj_node = get_node_safe(nodes, x + Stencil::cx(i), y + Stencil::cy(i));

        if (adj_node == to_u8(NodeId::FLUID))
        {
            nodes[idx] = to_u8(NodeId::DIRICHLET);
            break;
        }
    }
}

__global__ void init_valid_dirs(const uint8_t *__restrict__ nodes,
                                uint32_t *__restrict__ out_dirs)
{
    int x, y;
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    const uint8_t SOLID = to_u8(NodeId::SOLID);

    if (nodes[idx] == SOLID)
    {
        out_dirs[idx] = 0u;
        return;
    }

    uint32_t m = 0u;
    m |= (1u << 0);

#pragma unroll
    for (int i = 1; i < Stencil::Q; ++i)
    {
        const int cx = Stencil::cx(i);
        const int cy = Stencil::cy(i);

        const int xn = x + cx;
        const int yn = y + cy;

        if (get_node_safe(nodes, xn, yn) == SOLID)
            continue;

        if (get_node_safe(nodes, xn, y) == SOLID)
            continue;

        if (get_node_safe(nodes, x, yn) == SOLID)
            continue;

        m |= (1u << i);
    }

    out_dirs[idx] = m;
}

void build_tags(DomainTags &T)
{
    dim3 block(16, 16, 1);
    dim3 grid((NX + block.x - 1) / block.x,
              (NY + block.y - 1) / block.y, 1);

    annul_tags_kernel<<<grid, block>>>(T.d_node);
    CUDA_CHECK(cudaGetLastError());

    annul_tags_boundary_kernel<<<grid, block>>>(T.d_node);
    CUDA_CHECK(cudaGetLastError());

    init_valid_dirs<<<grid, block>>>(T.d_node, T.d_valid);
    CUDA_CHECK(cudaGetLastError());

    if (T.h_valid && T.h_node)
    {
        CUDA_CHECK(cudaMemcpy(T.h_valid, T.d_valid, T.bytes_valid, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(T.h_node, T.d_node, T.bytes_node, cudaMemcpyDeviceToHost));
    }
}
