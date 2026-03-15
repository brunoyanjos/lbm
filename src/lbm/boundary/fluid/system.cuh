// #pragma once

// #include "../../../core/types.cuh"
// #include "../common/incoming_moments.cuh"
// #include "../common/linear_row.cuh"
// #include "../common/equation_row.cuh"
// #include "layout.cuh"

// struct FluidSystem
// {
//     IncomingMoments<FluidIncomingList> incomings;

//     EquationRow<FluidVarList, FluidNonlinearList> rho;

//     LinearRow<FluidVarList> ux;
//     LinearRow<FluidVarList> uy;

//     LinearRow<FluidVarList> mxx;
//     LinearRow<FluidVarList> mxy;
//     LinearRow<FluidVarList> myy;

//     __device__ __forceinline__ void normalize_known()
//     {
//         const real_t inv_rho = r::one / incomings.template get<MomentId::rho>();

//         incomings.template get<MomentId::ux>() *= inv_rho;
//         incomings.template get<MomentId::uy>() *= inv_rho;
//         incomings.template get<MomentId::mxx>() *= inv_rho;
//         incomings.template get<MomentId::mxy>() *= inv_rho;
//         incomings.template get<MomentId::myy>() *= inv_rho;
//     }
// };