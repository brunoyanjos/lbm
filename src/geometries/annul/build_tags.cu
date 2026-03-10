#include "build_tags.cuh"

#include <cstdlib>

#include "properties.cuh"
#include "../../core/cuda_utils.cuh"
#include "../../core/indexing.cuh"
#include "../../core/math_utils.cuh"

#include "../../lbm/stencil_active.cuh"
#include "../../lbm/domain/mask_utils.cuh"

namespace ANNUL
{
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

            if (cx + cy == 4)
            {
                if (get_node_safe(nodes, x + cx / 2, y + cy / 2) == SOLID)
                    continue;
            }
            else if (cx + cy == 3)
            {
                if (get_node_safe(nodes, x + cx / 3, y + cy / 3) == SOLID || get_node_safe(nodes, x + cx * 2 / 3, y + cy * 2 / 3) == SOLID)
                    continue;
            }

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

        // annul_tags_boundary_kernel<<<grid, block>>>(T.d_node);
        // CUDA_CHECK(cudaGetLastError());

        init_valid_dirs<<<grid, block>>>(T.d_node, T.d_valid);
        CUDA_CHECK(cudaGetLastError());

        if (T.h_valid && T.h_node)
        {
            CUDA_CHECK(cudaMemcpy(T.h_valid, T.d_valid, T.bytes_valid, cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(T.h_node, T.d_node, T.bytes_node, cudaMemcpyDeviceToHost));
        }
    }
}