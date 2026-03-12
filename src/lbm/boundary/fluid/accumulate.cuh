#pragma once

#include "../../../core/types.cuh"
#include "../../../geometries/active_geometry.cuh"
#include "../../stencil_active.cuh"

#include "../common/moment_ids.cuh"
#include "../common/nonlinear_ids.cuh"

#include "system.cuh"
#include "../common/factors.cuh"

__device__ __forceinline__ void accumulate_incoming_fluid(
    FluidSystem &S,
    real_t pop_i,
    const Stencil::Basis &B,
    const Factors &F)
{
    S.incomings.template get<MomentId::rho>() += pop_i;
    S.incomings.template get<MomentId::ux>() += pop_i * B.cx;
    S.incomings.template get<MomentId::uy>() += pop_i * B.cy;

    S.incomings.template get<MomentId::mxx>() += pop_i * B.Hxx;
    S.incomings.template get<MomentId::mxy>() += pop_i * B.Hxy;
    S.incomings.template get<MomentId::myy>() += pop_i * B.Hyy;

    S.ux.template get<MomentId::rho>() += F.w * B.cx;
    S.ux.template get<MomentId::ux>() += F.Hx * B.cx;
    S.ux.template get<MomentId::uy>() += F.Hy * B.cx;
    S.ux.template get<MomentId::mxx>() += F.Hxx * B.cx;
    S.ux.template get<MomentId::mxy>() += F.Hxy * B.cx;
    S.ux.template get<MomentId::myy>() += F.Hyy * B.cx;

    S.uy.template get<MomentId::rho>() += F.w * B.cy;
    S.uy.template get<MomentId::ux>() += F.Hx * B.cy;
    S.uy.template get<MomentId::uy>() += F.Hy * B.cy;
    S.uy.template get<MomentId::mxx>() += F.Hxx * B.cy;
    S.uy.template get<MomentId::mxy>() += F.Hxy * B.cy;
    S.uy.template get<MomentId::myy>() += F.Hyy * B.cy;

    S.mxx.template get<MomentId::rho>() += F.w * B.Hxx;
    S.mxx.template get<MomentId::ux>() += F.Hx * B.Hxx;
    S.mxx.template get<MomentId::uy>() += F.Hy * B.Hxx;
    S.mxx.template get<MomentId::mxx>() += F.Hxx * B.Hxx;
    S.mxx.template get<MomentId::mxy>() += F.Hxy * B.Hxx;
    S.mxx.template get<MomentId::myy>() += F.Hyy * B.Hxx;

    S.mxy.template get<MomentId::rho>() += F.w * B.Hxy;
    S.mxy.template get<MomentId::ux>() += F.Hx * B.Hxy;
    S.mxy.template get<MomentId::uy>() += F.Hy * B.Hxy;
    S.mxy.template get<MomentId::mxx>() += F.Hxx * B.Hxy;
    S.mxy.template get<MomentId::mxy>() += F.Hxy * B.Hxy;
    S.mxy.template get<MomentId::myy>() += F.Hyy * B.Hxy;

    S.myy.template get<MomentId::rho>() += F.w * B.Hyy;
    S.myy.template get<MomentId::ux>() += F.Hx * B.Hyy;
    S.myy.template get<MomentId::uy>() += F.Hy * B.Hyy;
    S.myy.template get<MomentId::mxx>() += F.Hxx * B.Hyy;
    S.myy.template get<MomentId::mxy>() += F.Hxy * B.Hyy;
    S.myy.template get<MomentId::myy>() += F.Hyy * B.Hyy;
}

__device__ __forceinline__ void accumulate_outgoing_fluid(
    FluidSystem &S,
    const Factors &F)
{
    S.rho.lin.template get<MomentId::rho>() += F.w;
    S.rho.lin.template get<MomentId::ux>() += F.Hx;
    S.rho.lin.template get<MomentId::uy>() += F.Hy;

    S.rho.lin.template get<MomentId::mxx>() += (r::one - Geometry::OMEGA) * F.Hxx;
    S.rho.lin.template get<MomentId::mxy>() += (r::one - Geometry::OMEGA) * F.Hxy;
    S.rho.lin.template get<MomentId::myy>() += (r::one - Geometry::OMEGA) * F.Hyy;

    S.rho.nonlin.template get<NonlinearId::uxux>() += Geometry::OMEGA * F.Hxx;
    S.rho.nonlin.template get<NonlinearId::uxuy>() += Geometry::OMEGA * F.Hxy;
    S.rho.nonlin.template get<NonlinearId::uyuy>() += Geometry::OMEGA * F.Hyy;
}