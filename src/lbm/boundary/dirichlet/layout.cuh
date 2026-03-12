#pragma once

#include "../common/id_list.cuh"
#include "../common/moment_ids.cuh"
#include "../common/nonlinear_ids.cuh"

using DirichletVarList = IdList<
    MomentId::rho,
    MomentId::ux,
    MomentId::uy,
    MomentId::mxx,
    MomentId::mxy,
    MomentId::myy>;
// MomentId::mxxx,
// MomentId::mxxy,
// MomentId::mxyy,
// MomentId::myyy>;

using DirichletIncomingList = IdList<
    MomentId::rho,
    MomentId::mxx,
    MomentId::mxy,
    MomentId::myy>;
// MomentId::mxxx,
// MomentId::mxxy,
// MomentId::mxyy,
// MomentId::myyy>;

using DirichletEquationList = IdList<
    MomentId::mxx,
    MomentId::mxy,
    MomentId::myy>;
// MomentId::mxxx,
// MomentId::mxxy,
// MomentId::mxyy,
// MomentId::myyy>;

using DirichletNonlinearList = IdList<
    NonlinearId::uxux,
    NonlinearId::uxuy,
    NonlinearId::uyuy>;
// NonlinearId::uxuxux,
// NonlinearId::uxuxuy,
// NonlinearId::uxuyuy,
// NonlinearId::uyuyuy > ;