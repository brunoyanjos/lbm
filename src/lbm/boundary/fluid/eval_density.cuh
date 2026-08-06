// #pragma once

// #include "core/types.cuh"

// #include "lbm/boundary/fluid/accumulator.cuh"
// #include "lbm/moment/node_moments.cuh"

// namespace boundary::fluid
// {
//     template <int RegOrder, bool Rec, bool HighOrder>
//     [[nodiscard]] __device__ __forceinline__ real_t eval_density(const Accumulator<RegOrder, Rec> &acc,
//                                                                  const NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
//     {
//         real_t rho_denominator = acc.rho.rho +
//                                  M.ux * acc.rho.ux + M.uy * acc.rho.uy +
//                                  M.ux * M.ux * acc.rho.uxux + M.ux * M.uy * acc.rho.uxuy +
//                                  M.uy * M.uy * acc.rho.uyuy +
//                                  M.mxx * acc.rho.mxx + M.mxy * acc.rho.mxy +
//                                  M.myy * acc.rho.myy;

//         if constexpr (RegOrder >= 3)
//         {
//             rho_denominator += M.ux * M.ux * M.ux * acc.rho.uxuxux + M.ux * M.ux * M.uy * acc.rho.uxuxuy +
//                                M.ux * M.uy * M.uy * acc.rho.uxuyuy + M.uy * M.uy * M.uy * acc.rho.uyuyuy;

//             if constexpr (Rec)
//             {
//                 rho_denominator += M.ux * M.mxx * acc.rho.uxmxx + M.uy * M.mxx * acc.rho.uymxx +
//                                    M.ux * M.mxy * acc.rho.uxmxy + M.uy * M.mxy * acc.rho.uymxy +
//                                    M.ux * M.myy * acc.rho.uxmyy + M.uy * M.myy * acc.rho.uymyy;
//             }
//         }

//         const real_t inv_rho = r::one / rho_denominator;

//         return acc.in.rho * inv_rho;
//     }
// }