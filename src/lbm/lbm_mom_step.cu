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
    const uint32_t valid_ms = T.d_valid[idx];

    real_t pop[Stencil::Q];

    reconstruct_streamed_pop(pop, S, c, x, y);

    NodeMoments M{};

    if (node_id != to_u8(NodeId::FLUID))
    {
        bc_velocity(M, x, y);

        apply_boundary(pop, valid_ms, M);
    }
    else
    {
        if (is_full_mask(valid_ms))
        {
            evaluate_moments_from_pop(pop, M);
        }
        else
        {
            M.rho = S.d_rho[c][idxGlobal(x, y)] + RHO_0;
            M.ux = S.d_ux[c][idxGlobal(x, y)] * inv_scale_factor<MomentId::ux>();
            M.uy = S.d_uy[c][idxGlobal(x, y)] * inv_scale_factor<MomentId::uy>();
            M.mxx = S.d_mxx[c][idxGlobal(x, y)] * inv_scale_factor<MomentId::mxx>();
            M.mxy = S.d_mxy[c][idxGlobal(x, y)] * inv_scale_factor<MomentId::mxy>();
            M.myy = S.d_myy[c][idxGlobal(x, y)] * inv_scale_factor<MomentId::myy>();

            evaluate_fluid_node(pop, valid_ms, M);
        }
    }

    // 3) scale to the stored basis
    scale_to_stored_basis(M);

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
