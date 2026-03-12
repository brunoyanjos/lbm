#pragma once

#include "../../../core/types.cuh"
#include "../common/incoming_moments.cuh"
#include "../common/linear_row.cuh"
#include "../common/equation_row.cuh"
#include "layout.cuh"

struct FluidSystem2D
{
    IncomingMoments<FluidIncomingList2D> incomings;

    EquationRow<FluidVarList2D, FluidNonlinearList2D> rho;

    LinearRow<FluidVarList2D> ux;
    LinearRow<FluidVarList2D> uy;

    LinearRow<FluidVarList2D> mxx;
    LinearRow<FluidVarList2D> mxy;
    LinearRow<FluidVarList2D> myy;

    __device__ __forceinline__ void normalize_known()
    {
        const real_t inv_rho = r::one / incomings.template get<MomentId2D::rho>();

        incomings.template get<MomentId2D::ux>() *= inv_rho;
        incomings.template get<MomentId2D::uy>() *= inv_rho;
        incomings.template get<MomentId2D::mxx>() *= inv_rho;
        incomings.template get<MomentId2D::mxy>() *= inv_rho;
        incomings.template get<MomentId2D::myy>() *= inv_rho;
    }
};