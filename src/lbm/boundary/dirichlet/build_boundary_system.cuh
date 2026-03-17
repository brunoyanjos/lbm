#pragma once

#include <cstddef>

#include "core/types.cuh"

#include "lbm/boundary/common/system_data.cuh"
#include "lbm/boundary/common/accumulator.cuh"

template <std::size_t N, bool HasVelocity>
__device__ __forceinline__ void build_boundary_system(SystemData<N> &S, real_t ux, real_t uy, const MomentAccumulator<HasVelocity> &acc)
{
    const real_t rho_eq = acc.rho.rho + ux * acc.rho.ux + uy * acc.rho.uy +
                          ux * ux * acc.rho.uxux + ux * uy * acc.rho.uxuy + uy * uy * acc.rho.uyuy;

    // mxx
    S.coeff(0, 0) = acc.mxx.mxx - acc.rho.mxx * acc.in.mxx;
    S.coeff(0, 1) = acc.mxx.mxy - acc.rho.mxy * acc.in.mxx;
    S.coeff(0, 2) = acc.mxx.myy - acc.rho.myy * acc.in.mxx;
    S.b[0] = acc.in.mxx * rho_eq - (acc.mxx.rho + ux * acc.mxx.ux + uy * acc.mxx.uy);

    // mxy
    S.coeff(1, 0) = acc.mxy.mxx - acc.rho.mxx * acc.in.mxy;
    S.coeff(1, 1) = acc.mxy.mxy - acc.rho.mxy * acc.in.mxy;
    S.coeff(1, 2) = acc.mxy.myy - acc.rho.myy * acc.in.mxy;
    S.b[1] = acc.in.mxy * rho_eq - (acc.mxy.rho + ux * acc.mxy.ux + uy * acc.mxy.uy);

    // myy
    S.coeff(2, 0) = acc.myy.mxx - acc.rho.mxx * acc.in.myy;
    S.coeff(2, 1) = acc.myy.mxy - acc.rho.mxy * acc.in.myy;
    S.coeff(2, 2) = acc.myy.myy - acc.rho.myy * acc.in.myy;
    S.b[2] = acc.in.myy * rho_eq - (acc.myy.rho + ux * acc.myy.ux + uy * acc.myy.uy);
}