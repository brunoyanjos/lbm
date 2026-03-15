#pragma once

#include "../../core/types.cuh"
#include "../moment/moment_values.cuh"

template <auto Id, class MomentList>
__host__ __device__ __forceinline__
    real_t
    equilibrium_moment_value(const MomentValues<MomentList> &M)
{
    const real_t ux = M.template get<MomentId::ux>();
    const real_t uy = M.template get<MomentId::uy>();

    if constexpr (Id == MomentId::rho)
        return M.template get<MomentId::rho>();
    else if constexpr (Id == MomentId::ux)
        return ux;
    else if constexpr (Id == MomentId::uy)
        return uy;
    else if constexpr (Id == MomentId::mxx)
        return r::half * ux * ux;
    else if constexpr (Id == MomentId::mxy)
        return ux * uy;
    else if constexpr (Id == MomentId::myy)
        return r::half * uy * uy;
    else if constexpr (Id == MomentId::mxxy)
        return r::half * ux * ux * uy;
    else if constexpr (Id == MomentId::mxyy)
        return r::half * ux * uy * uy;
    else if constexpr (Id == MomentId::mxxx)
        return r::sixth * ux * ux * ux;
    else if constexpr (Id == MomentId::myyy)
        return r::sixth * uy * uy * uy;
    else
        static_assert(Id != Id, "MomentId sem equilíbrio associado.");
}