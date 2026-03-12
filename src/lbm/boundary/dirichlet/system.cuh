#pragma once

#include "../../../core/types.cuh"
#include "../common/incoming_moments.cuh"
#include "../common/linear_row.cuh"
#include "../common/equation_row.cuh"
#include "layout.cuh"

struct DirichletSystem2D
{
    IncomingMoments<DirichletIncomingList2D> incomings;

    EquationRow<DirichletVarList2D, DirichletNonlinearList2D> rho;

    LinearRow<DirichletVarList2D> mxx;
    LinearRow<DirichletVarList2D> mxy;
    LinearRow<DirichletVarList2D> myy;

    __device__ __forceinline__ void normalize_known()
    {
        const real_t inv_rho = r::one / incomings.template get<MomentId2D::rho>();

        incomings.template get<MomentId2D::mxx>() *= inv_rho;
        incomings.template get<MomentId2D::mxy>() *= inv_rho;
        incomings.template get<MomentId2D::myy>() *= inv_rho;
    }
};