// #pragma once

// #include "../../../core/types.cuh"
// #include "moment_ids.cuh"
// #include "factors.cuh"
// #include "../../stencil_active.cuh"

// __device__ __forceinline__ real_t moment_basis_value(MomentId id, const Stencil::Basis &B)
// {
//     switch (id)
//     {
//     case MomentId::rho:
//         return r::one;
//     case MomentId::ux:
//         return B.cx;
//     case MomentId::uy:
//         return B.cy;
//     case MomentId::mxx:
//         return B.Hxx;
//     case MomentId::mxy:
//         return B.Hxy;
//     case MomentId::myy:
//         return B.Hyy;

// #if LBM_HAS_REG3_AXIAL
//     case MomentId::mxxx:
//         return B.Hxxx;
//     case MomentId::myyy:
//         return B.Hyyy;
// #endif

// #if LBM_HAS_REG3_CROSS
//     case MomentId::mxxy:
//         return B.Hxxy;
//     case MomentId::mxyy:
//         return B.Hxyy;
// #endif

//     default:
//         return r::zero;
//     }
// }
