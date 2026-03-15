// #pragma once

// #include "../../../core/types.cuh"
// #include "../common/incoming_moments.cuh"
// #include "../common/linear_row.cuh"
// #include "../common/equation_row.cuh"
// #include "layout.cuh"

// struct DirichletSystem
// {
//     IncomingMoments<DirichletIncomingList> incomings;

//     EquationRow<DirichletVarList, DirichletNonlinearList> rho;

//     LinearRow<DirichletVarList> mxx;
//     LinearRow<DirichletVarList> mxy;
//     LinearRow<DirichletVarList> myy;

//     __device__ __forceinline__ void normalize_incomings()
//     {
//         const real_t inv_rho = r::one / incomings.template get<MomentId::rho>();

//         incomings.template get<MomentId::mxx>() *= inv_rho;
//         incomings.template get<MomentId::mxy>() *= inv_rho;
//         incomings.template get<MomentId::myy>() *= inv_rho;
//     }
// };