// #pragma once

// #include "../../../core/types.cuh"
// #include "system.cuh"
// #include "layout.cuh"
// #include "../common/linear_row.cuh"
// #include "../common/moment_ids.cuh"

// template <auto EqId>
// __device__ __forceinline__
//     LinearRow<DirichletVarList> &
//     dirichlet_row(DirichletSystem &S)
// {
//     if constexpr (EqId == MomentId::rho)
//         return S.rho;
//     else if constexpr (EqId == MomentId::mxx)
//         return S.mxx;
//     else if constexpr (EqId == MomentId::mxy)
//         return S.mxy;
//     else if constexpr (EqId == MomentId::myy)
//         return S.myy;
//     else
//         static_assert(EqId == MomentId::rho, "Unsupported EqId in dirichlet_row");
// }
