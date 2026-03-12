#pragma once

#include "../common/type_list.cuh"
#include "../common/moment_ids.cuh"
#include "../common/nonlinear_ids.cuh"

using DirichletVarList2D = IdList<
    MomentId2D::rho,
    MomentId2D::ux,
    MomentId2D::uy,
    MomentId2D::mxx,
    MomentId2D::mxy,
    MomentId2D::myy>;

using DirichletIncomingList2D = IdList<
    MomentId2D::rho,
    MomentId2D::mxx,
    MomentId2D::mxy,
    MomentId2D::myy>;

using DirichletNonlinearList2D = IdList<
    NonlinearId2D::uxux,
    NonlinearId2D::uxuy,
    NonlinearId2D::uyuy>;