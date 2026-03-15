// #pragma once

// #include "../../../core/types.cuh"
// #include "../common/moment_ids.cuh"
// #include "../common/nonlinear_ids.cuh"
// #include "system.cuh"

// __device__ __forceinline__ real_t recover_fluid_rho(
//     const FluidSystem &S,
//     real_t ux,
//     real_t uy,
//     real_t mxx,
//     real_t mxy,
//     real_t myy)
// {
//     const real_t denom =
//         S.rho.lin.template get<MomentId::rho>() +
//         ux * S.rho.lin.template get<MomentId::ux>() +
//         uy * S.rho.lin.template get<MomentId::uy>() +
//         mxx * S.rho.lin.template get<MomentId::mxx>() +
//         mxy * S.rho.lin.template get<MomentId::mxy>() +
//         myy * S.rho.lin.template get<MomentId::myy>() +
//         ux * ux * S.rho.nonlin.template get<NonlinearId::uxux>() +
//         ux * uy * S.rho.nonlin.template get<NonlinearId::uxuy>() +
//         uy * uy * S.rho.nonlin.template get<NonlinearId::uyuy>();

//     return S.incomings.template get<MomentId::rho>() / denom;
// }