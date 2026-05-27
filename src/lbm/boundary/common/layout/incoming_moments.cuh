#pragma once

#include "../../moment/moment_id.cuh"
#include "id_list.cuh"

template <bool HasVelocity>
struct LinearSystemListImpl;

template <>
struct LinearSystemListImpl<true>
{
    using type = IdList<
        MomentId::ux,
        MomentId::uy,
        MomentId::mxx,
        MomentId::mxy,
        MomentId::myy,
        >;
};

template <>
struct LinearSystemListImpl<false>
{
    using type = IdList<
        MomentId::mxx,
        MomentId::mxy,
        MomentId::myy,
        >;
};

template <bool HasVelocity>
using LinearSystemList = typename LinearSystemListImpl<HasVelocity>::type;
