#pragma once

#include "../meta/id_list.cuh"
#include "../meta/id_list_concat.cuh"
#include "../meta/id_list_if.cuh"
#include "../moment/moment_id.cuh"
#include "../regularization/traits/reg_traits.cuh"

using StateBaseFieldList = IdList<
    MomentId::rho,
    MomentId::ux,
    MomentId::uy,
    MomentId::mxx,
    MomentId::mxy,
    MomentId::myy>;

using StateReg3CrossFieldList = IdList<
    MomentId::mxxy,
    MomentId::mxyy>;

using StateReg3AxialFieldList = IdList<
    MomentId::mxxx,
    MomentId::myyy>;

using StateFieldListWithCross = IdListConcatT<
    StateBaseFieldList,
    IdListIfT<
        lbm_reg::Active::has_reg3_cross,
        StateReg3CrossFieldList,
        IdList<>>>;

using StateFieldList = IdListConcatT<
    StateFieldListWithCross,
    IdListIfT<
        lbm_reg::Active::has_reg3_axial,
        StateReg3AxialFieldList,
        IdList<>>>;