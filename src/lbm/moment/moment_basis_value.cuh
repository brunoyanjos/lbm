#pragma once

#include "../../core/types.cuh"
#include "../hermite/hermite.cuh"

template <auto Id>
__host__ __device__ __forceinline__
    real_t
    moment_basis_value(int icx, int icy)
{
    if constexpr (Id == MomentId::rho)
    {
        return r::one;
    }
    else if constexpr (Id == MomentId::ux)
    {
        return r_cast(icx);
    }
    else if constexpr (Id == MomentId::uy)
    {
        return r_cast(icy);
    }
    else
    {
        const real_t cx = r_cast(icx);
        const real_t cy = r_cast(icy);

        if constexpr (Id == MomentId::mxx)
            return Hermite::Hxx<Stencil>(cx);
        else if constexpr (Id == MomentId::mxy)
            return Hermite::Hxy(cx, cy);
        else if constexpr (Id == MomentId::myy)
            return Hermite::Hyy<Stencil>(cy);
        else if constexpr (Id == MomentId::mxxy)
            return Hermite::Hxxy<Stencil>(cx, cy);
        else if constexpr (Id == MomentId::mxyy)
            return Hermite::Hxyy<Stencil>(cx, cy);
        else if constexpr (Id == MomentId::mxxx)
            return Hermite::Hxxx<Stencil>(cx);
        else if constexpr (Id == MomentId::myyy)
            return Hermite::Hyyy<Stencil>(cy);
        else
            static_assert(Id != Id, "MomentId sem base associada em moment_basis_value.");
    }
}