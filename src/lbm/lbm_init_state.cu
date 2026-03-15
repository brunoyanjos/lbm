#include "lbm_init_state.cuh"

#include <cuda_runtime.h>

#include "../core/cuda_utils.cuh"
#include "../core/indexing.cuh"

#include "../geometries/active_geometry.cuh"

#include "equilibrium/lbm_equilibrium.cuh"
#include "stencil_active.cuh"

#include "layout/moment_eval.cuh"
#include "moment/moment_values.cuh"
#include "moment/evaluate_moments_from_pop.cuh"
#include "moment/scale_to_stored_basis.cuh"
#include "state/state_store.cuh"
#include "meta/for_each_id.cuh"

__global__ void init_on_device(LBMState S)
{
    int x, y;
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    const real_t rho = Geometry::RHO_0;
    const real_t ux = r::zero;
    const real_t uy = r::zero;

    real_t pop[Stencil::Q];
    equilibrium(pop, rho, ux, uy);

    MomentValues<MomentEvalList> M{};
    evaluate_moments_from_pop(pop, M);

    scale_to_stored_basis(M);

    const int c = S.cur;
    store_next_state(S, c, idx, M);
}

void init_state(LBMState &S, const CudaConfig &cfg)
{
    S.cur = 0;

    init_on_device<<<cfg.grid, cfg.block>>>(S);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}

struct HostUploader
{
    LBMState &S;
    int c;

    template <auto Id>
    __host__ void operator()()
    {
        CUDA_CHECK(cudaMemcpy(
            host_field<Id>(S),
            device_field<Id>(S)[c],
            S.bytes_field,
            cudaMemcpyDeviceToHost));
    }
};

void upload_state_to_host(LBMState &S)
{
    const int c = S.cur;
    for_each_id_host<StateFieldList>(HostUploader{S, c});
}