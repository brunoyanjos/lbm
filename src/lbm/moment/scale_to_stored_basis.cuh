#pragma once

#include "../../core/types.cuh"
#include "../stencil_active.cuh"
#include "../meta/id_list.cuh"
#include "../moment/moment_values.cuh"

template <auto Id>
__host__ __device__ __forceinline__
    real_t
    stored_basis_scale()
{
    if constexpr (Id == MomentId::rho)
    {
        return r::one;
    }
    else if constexpr (Id == MomentId::ux || Id == MomentId::uy)
    {
        return Stencil::as2;
    }
    else if constexpr (Id == MomentId::mxx || Id == MomentId::myy)
    {
        return r::half * Stencil::as4;
    }
    else if constexpr (Id == MomentId::mxy)
    {
        return Stencil::as4;
    }
    else if constexpr (Id == MomentId::mxxy || Id == MomentId::mxyy)
    {
        return r::half * Stencil::as6;
    }
    else if constexpr (Id == MomentId::mxxx || Id == MomentId::myyy)
    {
        return r::sixth * Stencil::as6;
    }
    else
    {
        static_assert(Id != Id, "MomentId sem fator associado em stored_basis_scale().");
    }
}

template <class MomentList>
struct StoredBasisScaler
{
    MomentValues<MomentList> &M;

    template <auto Id>
    __device__ __forceinline__ void operator()()
    {
        M.template get<Id>() *= stored_basis_scale<Id>();
    }
};

template <class MomentList>
__device__ __forceinline__ void scale_to_stored_basis(MomentValues<MomentList> &M)
{
    for_each_id<MomentList>(StoredBasisScaler<MomentList>{M});
}