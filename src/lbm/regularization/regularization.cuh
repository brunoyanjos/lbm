#pragma once

#include "../../core/types.cuh"
#include "../../core/indexing.cuh"

#include "../stencil_active.cuh"
#include "../state/lbm_state.cuh"
#include "../layout/pop_reconstruction_moments.cuh"
#include "accumulator.cuh"

#include "../meta/for_each_id.cuh"

__device__ __forceinline__ void reconstruct_streamed_pop(real_t *__restrict__ pop,
                                                         const LBMState &S,
                                                         int c,
                                                         int x,
                                                         int y)
{
#pragma unroll
        for (int i = 0; i < Stencil::Q; ++i)
        {
                const int icx = Stencil::cx(i);
                const int icy = Stencil::cy(i);

                const size_t n_idx = idxGlobalPeriodic(x - icx, y - icy);

                const real_t rho = device_field<MomentId::rho>(S)[c][n_idx] + Geometry::RHO_0;

                real_t f = r::one;

                for_each_id<PopReconMomentList>(PopReconstructionAccumulator{S, c, n_idx, icx, icy, f});

                pop[i] = Stencil::w(i) * rho * f;
        }
}