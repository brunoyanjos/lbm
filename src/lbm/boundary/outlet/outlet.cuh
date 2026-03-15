// #pragma once
// #include "../../../core/types.cuh"
// #include "../../../geometries/active_geometry.cuh"
// #include "../../../core/math_utils.cuh"
// #include "../../../core/linear_solver.cuh"
// #include "../../stencil_active.cuh"
// #include "../../domain/mask_utils.cuh"
// #include "../../state/lbm_state.cuh"
// #include <cstdint>
// #include <cstdio>

// __device__ __forceinline__ void apply_boundary_outlet(
//     real_t *__restrict__ pop,
//     uint32_t valid_mask,
//     real_t rho,
//     real_t ux, real_t uy,
//     real_t &mxx, real_t &mxy, real_t &myy)
// {
//     const uint32_t incoming_mask = mask_opp(valid_mask);

//     // soma das válidas
//     real_t rho_I = r::zero;

//     real_t mxx_I = r::zero;
//     real_t mxy_I = r::zero;
//     real_t myy_I = r::zero;

//     // from mxxI equation
//     real_t sum_Is_wi_Hxx = r::zero;

//     real_t sum_Is_wi_Hx_Hxx = r::zero;
//     real_t sum_Is_wi_Hy_Hxx = r::zero;

//     real_t sum_Is_wi_Hxx_Hxx = r::zero;
//     real_t sum_Is_wi_Hxy_Hxx = r::zero;
//     real_t sum_Is_wi_Hyy_Hxx = r::zero;

//     // from mxyI equation
//     real_t sum_Is_wi_Hxy = r::zero;

//     real_t sum_Is_wi_Hx_Hxy = r::zero;
//     real_t sum_Is_wi_Hy_Hxy = r::zero;

//     real_t sum_Is_wi_Hxx_Hxy = r::zero;
//     real_t sum_Is_wi_Hxy_Hxy = r::zero;
//     real_t sum_Is_wi_Hyy_Hxy = r::zero;

//     // from myyI equation
//     real_t sum_Is_wi_Hyy = r::zero;

//     real_t sum_Is_wi_Hx_Hyy = r::zero;
//     real_t sum_Is_wi_Hy_Hyy = r::zero;

//     real_t sum_Is_wi_Hxx_Hyy = r::zero;
//     real_t sum_Is_wi_Hxy_Hyy = r::zero;
//     real_t sum_Is_wi_Hyy_Hyy = r::zero;

// #pragma unroll
//     for (int i = 0; i < Stencil::Q; ++i)
//     {
//         real_t cx, cy, Hxx, Hxy, Hyy;
//         Stencil::basis2(i, cx, cy, Hxx, Hxy, Hyy);

//         const real_t Hx_factor = Stencil::w(i) * Stencil::as2 * cx;
//         const real_t Hy_factor = Stencil::w(i) * Stencil::as2 * cy;

//         const real_t Hxx_factor = Stencil::w(i) * Stencil::as4 * r::half * Hxx;
//         const real_t Hxy_factor = Stencil::w(i) * Stencil::as4 * r::half * Hxy;
//         const real_t Hyy_factor = Stencil::w(i) * Stencil::as4 * r::half * Hyy;

//         if (dir_valid(incoming_mask, i))
//         {
//             rho_I += pop[i];

//             mxx_I += pop[i] * Hxx;
//             mxy_I += pop[i] * Hxy;
//             myy_I += pop[i] * Hyy;

//             sum_Is_wi_Hxx += Stencil::w(i) * Hxx;
//             sum_Is_wi_Hx_Hxx += Hx_factor * Hxx;
//             sum_Is_wi_Hy_Hxx += Hy_factor * Hxx;
//             sum_Is_wi_Hxx_Hxx += Hxx_factor * Hxx;
//             sum_Is_wi_Hxy_Hxx += Hxy_factor * Hxx;
//             sum_Is_wi_Hyy_Hxx += Hyy_factor * Hxx;

//             sum_Is_wi_Hxy += Stencil::w(i) * Hxy;
//             sum_Is_wi_Hx_Hxy += Hx_factor * Hxy;
//             sum_Is_wi_Hy_Hxy += Hy_factor * Hxy;
//             sum_Is_wi_Hxx_Hxy += Hxx_factor * Hxy;
//             sum_Is_wi_Hxy_Hxy += Hxy_factor * Hxy;
//             sum_Is_wi_Hyy_Hxy += Hyy_factor * Hxy;

//             sum_Is_wi_Hyy += Stencil::w(i) * Hyy;
//             sum_Is_wi_Hx_Hyy += Hx_factor * Hyy;
//             sum_Is_wi_Hy_Hyy += Hy_factor * Hyy;
//             sum_Is_wi_Hxx_Hyy += Hxx_factor * Hyy;
//             sum_Is_wi_Hxy_Hyy += Hxy_factor * Hyy;
//             sum_Is_wi_Hyy_Hyy += Hyy_factor * Hyy;
//         }
//     }

//     const real_t inv_rho_I = real_t(1) / rho_I;

//     mxx_I *= inv_rho_I;
//     mxy_I *= inv_rho_I;
//     myy_I *= inv_rho_I;

//     constexpr int Nsys = 3;

//     real_t Aflat[Nsys * Nsys]{};

//     real_t b[Nsys]{};
//     real_t m[Nsys]{};

//     auto A = [&](int r, int c) -> real_t &
//     { return Aflat[r * Nsys + c]; };

//     const real_t inv_rho = r::one / rho;

//     // mxx
//     A(0, 0) = sum_Is_wi_Hxx_Hxx;
//     A(0, 1) = r::two * sum_Is_wi_Hxy_Hxx;
//     A(0, 2) = sum_Is_wi_Hyy_Hxx;
//     b[0] = rho_I * inv_rho * mxx_I - sum_Is_wi_Hxx -
//            ux * sum_Is_wi_Hx_Hxx - uy * sum_Is_wi_Hy_Hxx;

//     // mxy
//     A(1, 0) = sum_Is_wi_Hxx_Hxy;
//     A(1, 1) = r::two * sum_Is_wi_Hxy_Hxy;
//     A(1, 2) = sum_Is_wi_Hyy_Hxy;
//     b[1] = rho_I * inv_rho * mxy_I - sum_Is_wi_Hxy -
//            ux * sum_Is_wi_Hx_Hxy - uy * sum_Is_wi_Hy_Hxy;

//     // myy
//     A(2, 0) = sum_Is_wi_Hxx_Hyy;
//     A(2, 1) = r::two * sum_Is_wi_Hxy_Hyy;
//     A(2, 2) = sum_Is_wi_Hyy_Hyy;
//     b[2] = rho_I * inv_rho * myy_I - sum_Is_wi_Hyy -
//            ux * sum_Is_wi_Hx_Hyy - uy * sum_Is_wi_Hy_Hyy;

//     gaussianElimination<3>(Aflat, b, m);

//     mxx = m[0];
//     mxy = m[1];
//     myy = m[2];
// }