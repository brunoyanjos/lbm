#pragma once

#include "../../core/types.cuh"

#include "../stencil_active.cuh"
#include "../hermite/hermite.cuh"
#include "../layout/state_fields.cuh"
#include "../moment/moment_values.cuh"
#include "../moment/moment_basis_value.cuh"
#include "../meta/for_each_id.cuh"

#include "equilibrium_moment_value.cuh"

template <class MomentList>
struct EquilibriumAccumulator
{
    const MomentValues<MomentList> &M;
    int icx;
    int icy;
    real_t &f;

    template <auto Id>
    __device__ __forceinline__ void operator()()
    {
        const real_t meq = equilibrium_moment_value<Id>(M);
        f += meq * moment_basis_value<Id>(icx, icy);
    }
};

__device__ __forceinline__ void equilibrium(real_t *__restrict__ pop,
                                            real_t rho,
                                            real_t ux,
                                            real_t uy)
{
    MomentValues<StateFieldList> M{};

    M.template get<MomentId::ux>() = ux;
    M.template get<MomentId::uy>() = uy;
    M.template get<MomentId::mxx>() = r::zero;
    M.template get<MomentId::mxy>() = r::zero;
    M.template get<MomentId::myy>() = r::zero;

#pragma unroll
    for (int i = 0; i < Stencil::Q; ++i)
    {
        const int icx = Stencil::cx(i);
        const int icy = Stencil::cy(i);

        real_t f = r::one;

        for_each_id<StateFieldList>(EquilibriumAccumulator<StateFieldList>{M, icx, icy, f});

        pop[i] = Stencil::w(i) * rho * f;
    }
}