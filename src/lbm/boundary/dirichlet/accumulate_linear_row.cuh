#pragma once

#include "../../../core/types.cuh"
#include "../common/id_list.cuh"
#include "../common/moment_basis_value.cuh"
#include "../common/linear_factor_value.cuh"
#include "../common/linear_row.cuh"
#include "../common/factors.cuh"

template <class VarList>
struct OpMomLin
{
    LinearRow<VarList> &row;
    const Factors &F;
    real_t eq_basis;

    template <auto VarId>
    __device__ __forceinline__ void operator()()
    {
        row.template get<VarId>() += linear_factor_value(VarId, F) * eq_basis;
    }
};

template <auto EqId, class VarList>
__device__ __forceinline__ void accumulate_linear_row(
    LinearRow<VarList> &row,
    const Stencil::Basis &B,
    const Factors &F)
{
    const real_t eq_basis = moment_basis_value(EqId, B);

    for_each_id<VarList>(OpMomLin<VarList>{row, F, eq_basis});
}

template <class VarList>
struct OpRhoLin
{
    LinearRow<VarList> &row;
    const Factors &F;

    template <auto VarId>
    __device__ __forceinline__ void operator()()
    {
        real_t v = linear_factor_value(VarId, F);

        if constexpr (
            VarId == MomentId::rho ||
            VarId == MomentId::ux ||
            VarId == MomentId::uy)
        {
            row.template get<VarId>() += v;
        }
        else
        {
            row.template get<VarId>() += (r::one - Geometry::OMEGA) * v;
        }
    }
};

template <class VarList>
__device__ __forceinline__ void accumulate_rho_linear_row(
    LinearRow<VarList> &row,
    const Factors &F)
{
    for_each_id<VarList>(OpRhoLin<VarList>{row, F});
}