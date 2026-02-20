#include "lbm_mom_step.cuh"
#include "stencil_active.cuh"
#include "boundary/regularized_boundary_condition.cuh"
#include "collision/collision.cuh"
#include "moment_evaluation/moment_evaluation.cuh"
#include "moment_evaluation/moment_scaling.cuh"
#include "state/state_store.cuh"
#include "population/pop_reconstruction.cuh"

#include "../core/geometry.h"
#include "../core/physics.h"
#include "../core/indexing.cuh"
#include "../core/cuda_utils.cuh"

__global__ void lbm_mom_step_kernel(LBMState S, DomainTags T)
{
    int x, y;
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    const int c = S.cur;
    const int n = S.cur ^ 1;

    const uint8_t wall_id = T.d_wall[idx];
    const uint32_t valid_ms = T.d_valid[idx];

    real_t pop[Stencil::Q];

    // 1) reconstruct streamed populations from neighbor moments
    reconstruct_streamed_pop(pop, S, c, x, y);

    real_t rho, ux, uy, mxx, mxy, myy;

    evaluate_moments_from_pop(pop, rho, ux, uy, mxx, mxy, myy);

    if (wall_id != static_cast<uint8_t>(WallId::NONE))
    {
        ux = real_t(0);
        uy = real_t(0);

        if (wall_id == static_cast<uint8_t>(WallId::TOP))
            ux = U_LID;

        apply_boundary(pop, valid_ms, rho, ux, uy, mxx, mxy, myy, x, y, wall_id);
    }
    else
    {
        evaluate_moments_from_pop(pop, rho, ux, uy, mxx, mxy, myy);
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
