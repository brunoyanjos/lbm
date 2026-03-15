// #pragma once

// #include "../../../core/types.cuh"
// #include "../../../geometries/active_geometry.cuh"
// #include "../../stencil_active.cuh"
// #include "accumulate_incoming_moments.cuh"
// #include "accumulate_linear_row.cuh"
// #include "accumulate_nonlinear_row.cuh"
// #include "dirichlet_row.cuh"
// #include "system.cuh"
// #include "../common/factors.cuh"
// #include "../common/id_list.cuh"

// struct EqOp
// {
//     DirichletSystem &S;
//     const Stencil::Basis &B;
//     const Factors &F;

//     template <auto EqId>
//     __device__ __forceinline__ void operator()()
//     {
//         auto &row = dirichlet_row<EqId>(S);
//         accumulate_linear_row<EqId>(row, B, F);
//     }
// };

// __device__ __forceinline__ void accumulate_incoming_dirichlet(
//     DirichletSystem &S,
//     real_t pop_i,
//     const Stencil::Basis &B,
//     const Factors &F)
// {
//     accumulate_incoming_moments(S.incomings, pop_i, B);
//     for_each_id<DirichletEquationList>(EqOp{S, B, F});
// }

// __device__ __forceinline__ void accumulate_outgoing_dirichlet(
//     DirichletSystem &S,
//     const Factors &F)
// {
//     accumulate_rho_linear_row(S.rho.lin, F);
//     accumulate_rho_nonlinear_row(S.rho.nonlin, F);
// }