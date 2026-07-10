#pragma once

#include <cstddef>

#include "core/types.cuh"

#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/system_data.cuh"
#include "lbm/boundary/fluid/accumulator.cuh"

namespace boundary::fluid
{
    template <std::size_t N, std::size_t E>
    __device__ __forceinline__ void build_system(SystemData<N, E> &S,
                                                 const Accumulator &acc)
    {
        // uxI equation
        S.coeff(0, 0) = acc.ux.ux - acc.rho.ux * acc.in.ux;
        S.coeff(0, 1) = acc.ux.uy - acc.rho.uy * acc.in.ux;
        S.coeff(0, 2) = -acc.rho.uxux * acc.in.ux;
        S.coeff(0, 3) = -acc.rho.uxuy * acc.in.ux;
        S.coeff(0, 4) = -acc.rho.uyuy * acc.in.ux;
        S.coeff(0, 5) = acc.ux.mxx - acc.rho.mxx * acc.in.ux;
        S.coeff(0, 6) = acc.ux.mxy - acc.rho.mxy * acc.in.ux;
        S.coeff(0, 7) = acc.ux.myy - acc.rho.myy * acc.in.ux;

        S.b[0] = acc.rho.rho * acc.in.ux - acc.ux.rho;

        // uyI equation
        S.coeff(1, 0) = acc.uy.ux - acc.rho.ux * acc.in.uy;
        S.coeff(1, 1) = acc.uy.uy - acc.rho.uy * acc.in.uy;
        S.coeff(1, 2) = -acc.rho.uxux * acc.in.uy;
        S.coeff(1, 3) = -acc.rho.uxuy * acc.in.uy;
        S.coeff(1, 4) = -acc.rho.uyuy * acc.in.uy;
        S.coeff(1, 5) = acc.uy.mxx - acc.rho.mxx * acc.in.uy;
        S.coeff(1, 6) = acc.uy.mxy - acc.rho.mxy * acc.in.uy;
        S.coeff(1, 7) = acc.uy.myy - acc.rho.myy * acc.in.uy;

        S.b[1] = acc.rho.rho * acc.in.uy - acc.uy.rho;

        // mxxI equation
        S.coeff(2, 0) = acc.mxx.ux - acc.rho.ux * acc.in.mxx;
        S.coeff(2, 1) = acc.mxx.uy - acc.rho.uy * acc.in.mxx;
        S.coeff(2, 2) = -acc.rho.uxux * acc.in.mxx;
        S.coeff(2, 3) = -acc.rho.uxuy * acc.in.mxx;
        S.coeff(2, 4) = -acc.rho.uyuy * acc.in.mxx;
        S.coeff(2, 5) = acc.mxx.mxx - acc.rho.mxx * acc.in.mxx;
        S.coeff(2, 6) = acc.mxx.mxy - acc.rho.mxy * acc.in.mxx;
        S.coeff(2, 7) = acc.mxx.myy - acc.rho.myy * acc.in.mxx;

        S.b[2] = acc.rho.rho * acc.in.mxx - acc.mxx.rho;

        // mxyI equation
        S.coeff(3, 0) = acc.mxy.ux - acc.rho.ux * acc.in.mxy;
        S.coeff(3, 1) = acc.mxy.uy - acc.rho.uy * acc.in.mxy;
        S.coeff(3, 2) = -acc.rho.uxux * acc.in.mxy;
        S.coeff(3, 3) = -acc.rho.uxuy * acc.in.mxy;
        S.coeff(3, 4) = -acc.rho.uyuy * acc.in.mxy;
        S.coeff(3, 5) = acc.mxy.mxx - acc.rho.mxx * acc.in.mxy;
        S.coeff(3, 6) = acc.mxy.mxy - acc.rho.mxy * acc.in.mxy;
        S.coeff(3, 7) = acc.mxy.myy - acc.rho.myy * acc.in.mxy;

        S.b[3] = acc.rho.rho * acc.in.mxy - acc.mxy.rho;

        // myyI equation
        S.coeff(4, 0) = acc.myy.ux - acc.rho.ux * acc.in.myy;
        S.coeff(4, 1) = acc.myy.uy - acc.rho.uy * acc.in.myy;
        S.coeff(4, 2) = -acc.rho.uxux * acc.in.myy;
        S.coeff(4, 3) = -acc.rho.uxuy * acc.in.myy;
        S.coeff(4, 4) = -acc.rho.uyuy * acc.in.myy;
        S.coeff(4, 5) = acc.myy.mxx - acc.rho.mxx * acc.in.myy;
        S.coeff(4, 6) = acc.myy.mxy - acc.rho.mxy * acc.in.myy;
        S.coeff(4, 7) = acc.myy.myy - acc.rho.myy * acc.in.myy;

        S.b[4] = acc.rho.rho * acc.in.myy - acc.myy.rho;
    }
}
