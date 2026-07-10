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

    using FluidNonlinearMoments = IdList<
        NonlinearMomentId::uxux,
        NonlinearMomentId::uxuy,
        NonlinearMomentId::uyuy>;

}
