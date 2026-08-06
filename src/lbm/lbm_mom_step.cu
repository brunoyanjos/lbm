#include "core/geometry.h"
#include "core/physics.h"
#include "core/indexing.cuh"
#include "core/cuda_utils.cuh"

#include "lbm/lbm_mom_step.cuh"
#include "lbm/stencil_active.cuh"
#include "lbm/boundary/dirichlet/solver.cuh"
#include "lbm/boundary/fluid/solver.cuh"
#include "lbm/boundary/bc_velocity.cuh"
#include "lbm/collision/collision.cuh"
#include "lbm/interface/evaluate_normal.cuh"
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
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    const int c = S.cur;
    const int n = S.cur ^ 1;

    const uint8_t node_id = T.d_node[idx];
    const mask_t valid_ms = T.d_valid[idx];

    real_t popA[Stencil::Q];
    real_t popB[Stencil::Q];

    reconstruct_streamed_pop(popA, popB, S, c, x, y);

    NodeMoments M{};

    if (node_id != to_u8(NodeId::FLUID))
    {
        const mask_t incoming_mask = mask_opp(valid_ms);

        for (int i = 0; i < Stencil::Q; ++i)
        {
            if (!dir_valid(incoming_mask, i))
            {
                popA[i] = popA[Stencil::opp(i)];
                popB[i] = popB[Stencil::opp(i)];
            }
        }
    }

    evaluate_moments_from_pop(popA, popB, M);

    // 3) scale to the stored basis
    // scale_to_stored_basis(M);

    S.d_n[n][idx] = evaluate_normal(c, x, y, S);

    // 4) collide in moment space
    moment_space_collision(M);

    // 5) store next
    store_next_state(S, n, idx, M);
}

void lbm_mom_step(LBMState &S, const CudaConfig &cfg, const DomainTags &T)
{
    lbm_mom_step_kernel<<<cfg.grid, cfg.block>>>(S, T);
    CUDA_CHECK(cudaGetLastError());
}
