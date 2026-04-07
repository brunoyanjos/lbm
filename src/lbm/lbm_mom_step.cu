#include "core/geometry.h"
#include "core/physics.h"
#include "core/indexing.cuh"
#include "core/cuda_utils.cuh"

#include "lbm/lbm_mom_step.cuh"
#include "lbm/stencil_active.cuh"
#include "lbm/boundary/dirichlet/solver.cuh"
#include "lbm/boundary/bc_velocity.cuh"
#include "lbm/collision/collision.cuh"
#include "lbm/moment/moment_evaluation.cuh"
#include "lbm/moment/moment_scaling.cuh"
#include "lbm/moment/node_moments.cuh"
#include "lbm/moment/scale_factor.cuh"
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

    if (node_id == to_u8(NodeId::SOLID))
        return;

    real_t pop[Stencil::Q];

    reconstruct_streamed_pop(pop, S, c, x, y);

    NodeMoments M{};

    if (node_id != to_u8(NodeId::FLUID))
    {
        bc_velocity(x, y, M);

        apply_boundary(pop, valid_ms, M);
    }
    else
    {
        evaluate_moments_from_pop(pop, M);
    }

    scale_to_stored_basis(M);

    store_next_state(S, n, idx, M);
}

void lbm_mom_step(LBMState &S, const CudaConfig &cfg, const DomainTags &T)
{
    lbm_mom_step_kernel<<<cfg.grid, cfg.block>>>(S, T);
    CUDA_CHECK(cudaGetLastError());
}
