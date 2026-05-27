#pragma once

#include "lbm/boundary/common/id_list.cuh"

#include "lbm/moment/moment_id.cuh"

template <bool HasVelocity>
struct UnknownMomentListImpl;

template <>
struct UnknownMomentListImpl<true>
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
struct UnknownMomentListImpl<false>
{
    using type = IdList<
        MomentId::mxx,
        MomentId::mxy,
        MomentId::myy,
        >;
};

template <bool HasVelocity>
using UnknownMomentList = typename UnknownMomentListImpl<HasVelocity>::type;
