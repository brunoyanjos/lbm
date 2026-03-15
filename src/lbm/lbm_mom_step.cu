#include "lbm_mom_step.cuh"
#include "stencil_active.cuh"
#include "boundary/fluid_boundary_evaluation.cuh"
#include "collision/collision.cuh"
#include "moment/evaluate_moments_from_pop.cuh"
#include "moment/scale_to_stored_basis.cuh"
#include "state/state_store.cuh"
#include "regularization/regularization.cuh"

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

    real_t pop[Stencil::Q];

    reconstruct_streamed_pop(pop, S, c, x, y);

    MomentValues<MomentEvalList> M{};
    evaluate_moments_from_pop(pop, M);

    scale_to_stored_basis(M);
    moment_space_collision(M);
    store_next_state(S, n, idx, M);
}

void lbm_mom_step(LBMState &S, const CudaConfig &cfg, const DomainTags &T)
{
    lbm_mom_step_kernel<<<cfg.grid, cfg.block>>>(S, T);
    CUDA_CHECK(cudaGetLastError());
}
