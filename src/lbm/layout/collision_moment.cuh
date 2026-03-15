#pragma once

#include "../meta/id_list.cuh"
#include "../meta/id_list_concat.cuh"
#include "../meta/id_list_if.cuh"
#include "../moment/moment_id.cuh"
#include "../regularization/traits/reg_traits.cuh"

using CollisionMomentBaseList = IdList<
    MomentId::mxx,
    MomentId::mxy,
    MomentId::myy>;

using CollisionMomentReg3CrossList = IdList<
    MomentId::mxxy,
    MomentId::mxyy>;

using CollisionMomentReg3AxialList = IdList<
    MomentId::mxxx,
    MomentId::myyy>;

using CollisionMomentListWithCross = IdListConcatT<
    CollisionMomentBaseList,
    IdListIfT<
        lbm_reg::Active::has_reg3_cross,
        CollisionMomentReg3CrossList,
        IdList<>>>;

using CollisionMomentList = IdListConcatT<
    CollisionMomentListWithCross,
    IdListIfT<
        lbm_reg::Active::has_reg3_axial,
        CollisionMomentReg3AxialList,
        IdList<>>>;