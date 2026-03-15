#pragma once

#include "../meta/id_list.cuh"
#include "../meta/id_list_concat.cuh"
#include "../meta/id_list_if.cuh"

#include "../moment/moment_id.cuh"
#include "../regularization/traits/reg_traits.cuh"

using PopReconMomentBase = IdList<
    MomentId::ux,
    MomentId::uy,
    MomentId::mxx,
    MomentId::mxy,
    MomentId::myy>;

using PopReconMomentReg3Cross = IdList<
    MomentId::mxxy,
    MomentId::mxyy>;

using PopReconMomentReg3Axial = IdList<
    MomentId::mxxx,
    MomentId::myyy>;

using PopReconMomentListWithCross = IdListConcatT<
    PopReconMomentBase,
    IdListIfT<
        lbm_reg::Active::has_reg3_cross,
        PopReconMomentReg3Cross,
        IdList<>>>;

using PopReconMomentList = IdListConcatT<
    PopReconMomentListWithCross,
    IdListIfT<
        lbm_reg::Active::has_reg3_axial,
        PopReconMomentReg3Axial,
        IdList<>>>;