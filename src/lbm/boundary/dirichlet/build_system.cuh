// #pragma once

// #include "../../../core/types.cuh"
// #include "../../../core/linear_solver.cuh"
// #include "system.cuh"

// template <int N>
// __device__ __forceinline__ void build_dirichlet_system(
//     const DirichletSystem &S,
//     real_t ux,
//     real_t uy,
//     LinearSystem<N> &L)
// {
//     const real_t rho_part =
//         S.rho.lin.template get<MomentId::rho>() +
//         ux * S.rho.lin.template get<MomentId::ux>() +
//         uy * S.rho.lin.template get<MomentId::uy>() +
//         ux * ux * S.rho.nonlin.template get<NonlinearId::uxux>() +
//         ux * uy * S.rho.nonlin.template get<NonlinearId::uxuy>() +
//         uy * uy * S.rho.nonlin.template get<NonlinearId::uyuy>();

//     const real_t mxxI = S.incomings.template get<MomentId::mxx>();
//     const real_t mxyI = S.incomings.template get<MomentId::mxy>();
//     const real_t myyI = S.incomings.template get<MomentId::myy>();

//     L.coeff(0, 0) = S.mxx.template get<MomentId::mxx>() - S.rho.lin.template get<MomentId::mxx>() * mxxI;
//     L.coeff(0, 1) = S.mxx.template get<MomentId::mxy>() - S.rho.lin.template get<MomentId::mxy>() * mxxI;
//     L.coeff(0, 2) = S.mxx.template get<MomentId::myy>() - S.rho.lin.template get<MomentId::myy>() * mxxI;
//     L.b[0] = mxxI * rho_part - (S.mxx.template get<MomentId::rho>() +
//                                 ux * S.mxx.template get<MomentId::ux>() + uy * S.mxx.template get<MomentId::uy>());

//     L.coeff(1, 0) = S.mxy.template get<MomentId::mxx>() - S.rho.lin.template get<MomentId::mxx>() * mxyI;
//     L.coeff(1, 1) = S.mxy.template get<MomentId::mxy>() - S.rho.lin.template get<MomentId::mxy>() * mxyI;
//     L.coeff(1, 2) = S.mxy.template get<MomentId::myy>() - S.rho.lin.template get<MomentId::myy>() * mxyI;
//     L.b[1] = mxyI * rho_part - (S.mxy.template get<MomentId::rho>() +
//                                 ux * S.mxy.template get<MomentId::ux>() + uy * S.mxy.template get<MomentId::uy>());

//     L.coeff(2, 0) = S.myy.template get<MomentId::mxx>() - S.rho.lin.template get<MomentId::mxx>() * myyI;
//     L.coeff(2, 1) = S.myy.template get<MomentId::mxy>() - S.rho.lin.template get<MomentId::mxy>() * myyI;
//     L.coeff(2, 2) = S.myy.template get<MomentId::myy>() - S.rho.lin.template get<MomentId::myy>() * myyI;
//     L.b[2] = myyI * rho_part - (S.myy.template get<MomentId::rho>() +
//                                 ux * S.myy.template get<MomentId::ux>() + uy * S.myy.template get<MomentId::uy>());
// }