#pragma once

#include "lbm/boundary/common/id_list.cuh"

#include "lbm/moment/moment_id.cuh"

namespace boundary::fluid
{

    using UnknownMomentList = IdList<
        MomentId::ux,
        MomentId::uy,
        MomentId::mxx,
        MomentId::mxy,
        MomentId::myy,
        >;

    using SecondOrderNonLinearMomentList = IdList<
        NonlinearMomentId::uxux,
        NonlinearMomentId::uxuy,
        NonlinearMomentId::uyuy>;

    using ThirdOrderNonLinearMomentList = IdList<
        NonlinearMomentId::uxux,
        NonlinearMomentId::uxuy,
        NonlinearMomentId::uyuy,
        NonlinearMomentId::uxuxux,
        NonlinearMomentId::uxuxuy,
        NonlinearMomentId::uxuyuy,
        NonlinearMomentId::uyuyuy>;

    using ThirdOrderRecNonLinearMomentList = IdList<
        NonlinearMomentId::uxux,
        NonlinearMomentId::uxuy,
        NonlinearMomentId::uyuy,
        NonlinearMomentId::uxuxux,
        NonlinearMomentId::uxuxuy,
        NonlinearMomentId::uxuyuy,
        NonlinearMomentId::uyuyuy,
        NonlinearMomentId::uxmxx,
        NonlinearMomentId::uymxx,
        NonlinearMomentId::uxmxy,
        NonlinearMomentId::uymxy,
        NonlinearMomentId::uxmyy,
        NonlinearMomentId::uymyy>;

    template <int RegOrder, bool Rec>
    struct FluidNonlinearMomentsFor
    {
        static_assert(RegOrder == 2 || RegOrder == 3,
                      "Unsupported fluid boundary regularization order");

        using type = SecondOrderNonLinearMomentList;
    };

    template <>
    struct FluidNonlinearMomentsFor<3, false>
    {
        using type = ThirdOrderNonLinearMomentList;
    };

    template <>
    struct FluidNonlinearMomentsFor<3, true>
    {
        using type = ThirdOrderRecNonLinearMomentList;
    };

    template <int RegOrder, bool Rec>
    using FluidNonlinearMoments =
        typename FluidNonlinearMomentsFor<RegOrder, Rec>::type;

}
