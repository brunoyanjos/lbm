#include "core/geometry.h"
#include "core/physics.h"
#include "core/indexing.cuh"
#include "core/cuda_utils.cuh"
#include "app/cuda_config.cuh"

#include "lbm/lbm_mom_step.cuh"
#include "lbm/stencil_active.cuh"
#include "lbm/boundary/dirichlet/solver.cuh"
#include "lbm/boundary/fluid/solver.cuh"
#include "lbm/boundary/bc_velocity.cuh"
#include "lbm/collision/collision.cuh"
#include "lbm/moment/moment_evaluation.cuh"
#include "lbm/moment/moment_scaling.cuh"
#include "lbm/moment/node_moments.cuh"
#include "lbm/moment/scale_factor.cuh"
#include "lbm/state/state_load.cuh"
#include "lbm/state/state_store.cuh"
#include "lbm/population/pop_reconstruction.cuh"

__global__ void lbm_mom_step_kernel(LBMState S, DomainTags T)
{
    int x, y;
    const size_t idx = idxThreadLocalInterior2D(x, y, S.domain);
    if (idx == INVALID_INDEX)
        return;

    const int c = S.cur;
    const int n = S.cur ^ 1;

    const uint8_t node_id = T.d_node[idx];
    const mask_t valid_ms = T.d_valid[idx];

    real_t pop[Stencil::Q];

    reconstruct_streamed_pop(pop, S, c, x, y);

    NodeMoments M{};

    if (node_id != to_u8(NodeId::FLUID))
    {
        bc_velocity(M, x, y);

        boundary::dirichlet::apply_boundary(pop, valid_ms, M);
    }
    else
    {
        if (is_full_mask(valid_ms))
        {
            evaluate_moments_from_pop(pop, M);
        }
        else
        {
            load_state_moments(S, c, idx, M);

            boundary::fluid::apply_boundary(pop, valid_ms, M);
        }
    }

    // 3) scale to the stored basis
    scale_to_stored_basis(M);

    // 4) collide in moment space
    moment_space_collision(M);

    // 5) store next
    store_next_state(S, n, idx, M);
}

void lbm_mom_step(std::vector<LBMState> &S,
                  const std::vector<DomainTags> &T,
                  const app::RunContext &ctx)
{
    for (size_t i = 0; i < S.size(); ++i)
    {
        CUDA_CHECK(cudaSetDevice(ctx.partitions[i].device_id));
        const CudaConfig cfg = make_config(S[i].domain.nx, S[i].domain.local_ny);

        lbm_mom_step_kernel<<<cfg.grid, cfg.block>>>(S[i], T[i]);
        CUDA_CHECK(cudaGetLastError());
    }

    for (const auto &partition : ctx.partitions)
    {
        CUDA_CHECK(cudaSetDevice(partition.device_id));
        CUDA_CHECK(cudaDeviceSynchronize());
    }

    for (LBMState &state : S)
        state.cur ^= 1;
}
