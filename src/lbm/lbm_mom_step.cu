#include "lbm_mom_step.cuh"
#include "stencil_active.cuh"
#include "boundary/dirichlet/solver.cuh"
#include "boundary/fluid/solver.cuh"
#include "boundary/fluid_boundary_evaluation.cuh"
#include "collision/collision.cuh"
#include "moment_evaluation/moment_evaluation.cuh"
#include "moment_evaluation/moment_scaling.cuh"
#include "state/state_store.cuh"
#include "population/pop_reconstruction.cuh"

#include "../geometries/active_geometry.cuh"

#include "../core/indexing.cuh"
#include "../core/cuda_utils.cuh"
#include "../core/cuda_utils.cuh"

__global__ void lbm_mom_step_kernel(LBMState S, DomainTags T)
{
    int x, y;
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    const int c = S.cur;
    const int n = S.cur ^ 1;

    const uint8_t wall_id = T.d_node[idx];
    const uint32_t valid_ms = T.d_valid[idx];

    if (wall_id == to_u8(NodeId::SOLID))
        return;

    real_t pop[Stencil::Q];

    // 1) reconstruct streamed populations from neighbor moments
    reconstruct_streamed_pop(pop, S, c, x, y);

    real_t rho, ux, uy, mxx, mxy, myy;

    if (wall_id != static_cast<uint8_t>(NodeId::FLUID))
    {
        Geometry::bc_velocity(x, y, ux, uy);

        apply_boundary_dirichlet(pop, valid_ms, rho, ux, uy, mxx, mxy, myy);
    }
    else
    {
        if (is_full_mask(valid_ms))
        {
            evaluate_moments_from_pop(pop, rho, ux, uy, mxx, mxy, myy);
        }
        else
        {
            rho = S.d_rho[c][idx] + Geometry::RHO_0;
            ux = S.d_ux[c][idx] / Stencil::as2;
            uy = S.d_uy[c][idx] / Stencil::as2;
            mxx = S.d2.mxx[c][idx] / (Stencil::as4 * r::half);
            mxy = S.d2.mxy[c][idx] / Stencil::as4;
            myy = S.d2.myy[c][idx] / (Stencil::as4 * r::half);

            evaluate_fluid_node(pop, valid_ms, rho, ux, uy, mxx, mxy, myy);
        }
    }

    // 3) scale to the stored basis
    scale_to_stored_basis(ux, uy, mxx, mxy, myy);

    // 4) collide in moment space
    moment_space_collision(ux, uy, mxx, mxy, myy);

    // 5) store next
    store_next_state(S, n, idx, rho, ux, uy, mxx, mxy, myy);
}

void lbm_mom_step(LBMState &S, const CudaConfig &cfg, const DomainTags &T)
{
    lbm_mom_step_kernel<<<cfg.grid, cfg.block>>>(S, T);
    CUDA_CHECK(cudaGetLastError());
}
