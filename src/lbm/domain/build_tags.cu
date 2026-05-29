#include "build_tags.cuh"

#include "app/cuda_config.cuh"

#include "core/geometry.h"
#include "core/indexing.cuh"
#include "core/cuda_utils.cuh"

#include "lbm/stencil_active.cuh"
#include "lbm/domain/mask_utils.cuh"

#include <cstdlib>
#include <new>
#include <vector>

DomainTags domain_tags_allocate(const LocalDomain &domain)
{
    DomainTags T{};
    T.domain = domain;
    T.N = domain.N;
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
                                          uint8_t *__restrict__ node,
                                          LocalDomain domain)
{
    int x, y;
    const size_t idx = idxThreadLocalInterior2D(x, y, domain);
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
    CudaConfig cfg = make_config(T.domain.nx, T.domain.local_ny);

    cavity_square_tags_kernel<<<cfg.grid, cfg.block>>>(T.d_valid, T.d_node, T.domain);
    CUDA_CHECK(cudaGetLastError());

    if (T.h_valid && T.h_node)
    {
        CUDA_CHECK(cudaMemcpy(T.h_valid, T.d_valid, T.bytes_valid, cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpy(T.h_node, T.d_node, T.bytes_node, cudaMemcpyDeviceToHost));
    }
}

[[nodiscard]] __host__ std::vector<DomainTags> allocate_partition_tags(const app::RunContext &ctx)
{
    std::vector<DomainTags> tags;
    tags.reserve(ctx.partitions.size());

    for (const auto &partition : ctx.partitions)
    {
        CUDA_CHECK(cudaSetDevice(partition.device_id));
        const LocalDomain domain = make_local_domain(partition.y_begin,
                                                     partition.y_end,
                                                     partition.halo);

        tags.push_back(domain_tags_allocate(domain));
        build_tags(tags.back());
    }

    return tags;
}

__host__ void free_partition_tags(std::vector<DomainTags> &tags, const app::RunContext &ctx)
{
    for (size_t i = 0; i < tags.size(); ++i)
    {
        CUDA_CHECK(cudaSetDevice(ctx.partitions[i].device_id));
        domain_tags_free(tags[i]);
    }

    tags.clear();
}
