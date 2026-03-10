#include "build_tags.cuh"

#include "properties.cuh"
#include "../../core/indexing.cuh"
#include "../../lbm/stencil_active.cuh"
#include "../../core/cuda_utils.cuh"

#include <cstdlib>
#include <new>

namespace COUETTE
{
    __global__ void couette_tags_kernel(uint32_t *__restrict__ valid,
                                        uint8_t *__restrict__ wall)
    {
        int x, y;
        const size_t idx = idxThreadGlobal2D(x, y);
        if (idx == INVALID_INDEX)
            return;

        const bool on_bottom = (y == 0);
        const bool on_top = (y == NY - 1);

        uint8_t wid = static_cast<uint8_t>(NodeId::FLUID);
        if (on_bottom || on_top)
            wid = static_cast<uint8_t>(NodeId::DIRICHLET);

        wall[idx] = wid;

        uint32_t m = 0u;
        m |= (1u << 0);

#pragma unroll
        for (int i = 1; i < Stencil::Q; ++i)
        {
            int xn = x + Stencil::cx(i);
            int yn = y + Stencil::cy(i);

            // periódico em x
            if (xn < 0)
                xn += NX;
            else if (xn >= NX)
                xn -= NX;

            // em y, corta fora do domínio
            if (yn < 0 || yn >= NY)
                continue;

            m |= (1u << i);
        }

        valid[idx] = m;
    }

    void build_tags(DomainTags &T)
    {
        dim3 block(16, 16, 1);
        dim3 grid((NX + block.x - 1) / block.x,
                  (NY + block.y - 1) / block.y, 1);

        couette_tags_kernel<<<grid, block>>>(T.d_valid, T.d_node);
        CUDA_CHECK(cudaGetLastError());

        if (T.h_valid && T.h_node)
        {
            CUDA_CHECK(cudaMemcpy(T.h_valid, T.d_valid, T.bytes_valid, cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(T.h_node, T.d_node, T.bytes_node, cudaMemcpyDeviceToHost));
        }
    }
}
