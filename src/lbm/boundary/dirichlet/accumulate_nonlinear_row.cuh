#pragma once

#include "../../../core/types.cuh"
#include "../common/id_list.cuh"
#include "../common/moment_basis_value.cuh"
#include "../common/nonlinear_factor_value.cuh"
#include "../common/nonlinear_row.cuh"
#include "../common/factors.cuh"

template <class NonlinearList>
struct OpNonLinear
{
    NonlinearRow<NonlinearList> &row;
    const Factors &F;

    template <auto Id>
    __device__ __forceinline__ void operator()()
    {
        row.template get<Id>() += Geometry::OMEGA * nonlinear_factor_value(Id, F);
    }
};

template <class NonlinearList>
__device__ __forceinline__ void accumulate_rho_nonlinear_row(
    NonlinearRow<NonlinearList> &row,
    const Factors &F)
{
    for_each_id<NonlinearList>(OpNonLinear<NonlinearList>{row, F});
}