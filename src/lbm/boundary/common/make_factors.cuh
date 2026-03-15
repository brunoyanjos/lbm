// #pragma once

// #include "../../../core/types.cuh"
// #include "../../stencil_active.cuh"
// #include "factors.cuh"

// __device__ __forceinline__ Factors
// make_factors(const Stencil::Basis &B, real_t w)
// {
//     Factors F{};

//     F.w = w;
//     F.Hx = w * Stencil::as2 * B.cx;
//     F.Hy = w * Stencil::as2 * B.cy;
//     F.Hxx = w * Stencil::as4 * r::half * B.Hxx;
//     F.Hxy = w * Stencil::as4 * B.Hxy;
//     F.Hyy = w * Stencil::as4 * r::half * B.Hyy;

// #if LBM_HAS_REG3_CROSS
//     F.Hxxy = w * Stencil::as6 * r::half * B.Hxxy;
//     F.Hxyy = w * Stencil::as6 * r::half * B.Hxyy;
// #endif

// #if LBM_HAS_REG3_AXIAL
//     F.Hxxx = w * Stencil::as6 * r::sixth * B.Hxxx;
//     F.Hyyy = w * Stencil::as6 * r::sixth * B.Hyyy;
// #endif

//     return F;
// }