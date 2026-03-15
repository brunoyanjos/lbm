// #pragma once

// #include "../../../core/lbm_features.cuh"
// #include "../common/id_list.cuh"
// #include "../common/moment_ids.cuh"
// #include "../common/nonlinear_ids.cuh"

// using DirichletVarList = IdList<
//     MomentId::rho,
//     MomentId::ux,
//     MomentId::uy,
//     MomentId::mxx,
//     MomentId::mxy,
//     MomentId::myy

// #if LBM_HAS_REG3_CROSS
//     ,
//     MomentId::mxxx, MomentId::mxxy
// #endif

// #if LBM_HAS_REG3_AXIAL
//     ,
//     MomentId::mxyy, MomentId::myyy
// #endif
//     >;

// using DirichletIncomingList = IdList<
//     MomentId::rho,
//     MomentId::mxx,
//     MomentId::mxy,
//     MomentId::myy
// #if LBM_HAS_REG3_CROSS
//     ,
//     MomentId::mxxy, MomentId::mxyy
// #endif

// #if LBM_HAS_REG3_AXIAL
//     ,
//     MomentId::mxxx, MomentId::myyy
// #endif
//     >;

// using DirichletEquationList = IdList<
//     MomentId::mxx,
//     MomentId::mxy,
//     MomentId::myy,
// #if LBM_HAS_REG3_CROSS
//     ,
//     MomentId::mxxy, MomentId::mxyy
// #endif

// #if LBM_HAS_REG3_AXIAL
//     ,
//     MomentId::mxxx, MomentId::myyy
// #endif
//     >;

// using DirichletNonlinearList = IdList<
//     NonlinearId::uxux,
//     NonlinearId::uxuy,
//     NonlinearId::uyuy
// #if LBM_HAS_REG3_CROSS
//     ,
//     NonlinearId::uxuxuy, NonlinearId::uxuyuy
// #endif

// #if LBM_HAS_REG3_AXIAL
//     ,
//     NonlinearId::uxuxux, NonlinearId::uyuyuy
// #endif
//     >;