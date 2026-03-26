#pragma once

#include <cstddef>

#include "core/types.cuh"

#include "lbm/boundary/common/system_data.cuh"
#include "lbm/boundary/dirichlet/dirichlet_accumulator.cuh"

template <std::size_t N>
__device__ __forceinline__ void build_dirichlet_system(SystemData<N> &S, const DirichletAccumulator &acc)
{
    S.coeff(0, 0) = acc.mxx.mxx - acc.rho.mxx * acc.in.mxx;
    S.coeff(0, 1) = acc.mxx.mxy - acc.rho.mxy * acc.in.mxx;
    S.coeff(0, 2) = acc.mxx.myy - acc.rho.myy * acc.in.mxx;

    S.b[0] = acc.rho.constant * acc.in.mxx - acc.mxx.constant;

    S.coeff(1, 0) = acc.mxy.mxx - acc.rho.mxx * acc.in.mxy;
    S.coeff(1, 1) = acc.mxy.mxy - acc.rho.mxy * acc.in.mxy;
    S.coeff(1, 2) = acc.mxy.myy - acc.rho.myy * acc.in.mxy;

    S.b[1] = acc.rho.constant * acc.in.mxy - acc.mxy.constant;

    S.coeff(2, 0) = acc.myy.mxx - acc.rho.mxx * acc.in.myy;
    S.coeff(2, 1) = acc.myy.mxy - acc.rho.mxy * acc.in.myy;
    S.coeff(2, 2) = acc.myy.myy - acc.rho.myy * acc.in.myy;

    S.b[2] = acc.rho.constant * acc.in.myy - acc.myy.constant;
}