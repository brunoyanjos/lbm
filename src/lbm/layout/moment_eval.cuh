#pragma once

#include "../meta/id_list.cuh"
#include "../meta/id_list_concat.cuh"
#include "../meta/id_list_if.cuh"
#include "../moment/moment_id.cuh"
#include "../regularization/traits/reg_traits.cuh"

using MomentEvalBaseList = IdList<
    MomentId::rho,
    MomentId::ux,
    MomentId::uy,
    MomentId::mxx,
    MomentId::mxy,
    MomentId::myy>;

using MomentEvalReg3Cross = IdList<
    MomentId::mxxy,
    MomentId::mxyy>;

using MomentEvalReg3Axial = IdList<
    MomentId::mxxx,
    MomentId::myyy>;

using MomentEvalListWithCross = IdListConcatT<
    MomentEvalBaseList,
    IdListIfT<
        lbm_reg::Active::has_reg3_cross,
        MomentEvalReg3Cross,
        IdList<>>>;

using MomentEvalList = IdListConcatT<
    MomentEvalListWithCross,
    IdListIfT<
        lbm_reg::Active::has_reg3_axial,
        MomentEvalReg3Axial,
        IdList<>>>;