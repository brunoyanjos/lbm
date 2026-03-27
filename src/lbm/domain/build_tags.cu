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

    const real_t dx = r_cast(x) - x_center;
    const real_t dy = r_cast(y) - y_center;

    const real_t radius = r_sqrt(dx * dx + dy * dy);

    nodes[idx] = to_u8(NodeId::FLUID);

    if (radius < R_IN || radius > R_OUT)
        nodes[idx] = to_u8(NodeId::SOLID);
}

__device__ __forceinline__ uint8_t boundary_code_from_neighbors(
    uint8_t node_1, uint8_t node_2, uint8_t node_3, uint8_t node_4,
    uint8_t node_5, uint8_t node_6, uint8_t node_7, uint8_t node_8)
{
    const uint8_t FLUID = to_u8(NodeId::FLUID);

    const uint8_t bit_1 = (node_3 == FLUID || node_4 == FLUID || node_7 == FLUID) ? 1u : 0u;
    const uint8_t bit_2 = (node_1 == FLUID || node_4 == FLUID || node_8 == FLUID) ? 1u : 0u;
    const uint8_t bit_4 = (node_2 == FLUID || node_3 == FLUID || node_6 == FLUID) ? 1u : 0u;
    const uint8_t bit_8 = (node_1 == FLUID || node_2 == FLUID || node_5 == FLUID) ? 1u : 0u;

    return static_cast<uint8_t>(bit_1 | (bit_2 << 1) | (bit_4 << 2) | (bit_8 << 3));
}

__global__ void annul_boundary_id_kernel(uint8_t *__restrict__ nodes)
{
    int x, y;
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    if (nodes[idx] != to_u8(NodeId::SOLID))
        return;

    const uint8_t node_1 = get_node_safe(nodes, x + 1, y);
    const uint8_t node_2 = get_node_safe(nodes, x, y + 1);
    const uint8_t node_3 = get_node_safe(nodes, x - 1, y);
    const uint8_t node_4 = get_node_safe(nodes, x, y - 1);
    const uint8_t node_5 = get_node_safe(nodes, x + 1, y + 1);
    const uint8_t node_6 = get_node_safe(nodes, x - 1, y + 1);
    const uint8_t node_7 = get_node_safe(nodes, x - 1, y - 1);
    const uint8_t node_8 = get_node_safe(nodes, x + 1, y - 1);

    const uint8_t bc = boundary_code_from_neighbors(
        node_1, node_2, node_3, node_4,
        node_5, node_6, node_7, node_8);

    if (bc != 0u)
        nodes[idx] = bc;
}

void build_tags(DomainTags &T)
{
    dim3 block(16, 16, 1);
    dim3 grid((NX + block.x - 1) / block.x,
              (NY + block.y - 1) / block.y, 1);

    annul_tags_kernel<<<grid, block>>>(T.d_node);
    CUDA_CHECK(cudaGetLastError());

    annul_boundary_id_kernel<<<grid, block>>>(T.d_node);
    CUDA_CHECK(cudaGetLastError());

    if (T.h_node)
        CUDA_CHECK(cudaMemcpy(T.h_node, T.d_node, T.bytes_node, cudaMemcpyDeviceToHost));
}
