#pragma once

#include <cstddef>

#include "core/types.cuh"

#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/system_data.cuh"
#include "lbm/boundary/fluid/accumulator.cuh"

namespace boundary::fluid
{
    template <std::size_t N, std::size_t E, int RegOrder, bool Rec>
    __device__ __forceinline__ void build_system(SystemData<N, E> &S,
                                                 const Accumulator<RegOrder, Rec> &acc)
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

        if constexpr (RegOrder >= 3)
        {
            S.coeff(0, 8) = acc.ux.uxuxux - acc.rho.uxuxux * acc.in.ux;
            S.coeff(0, 9) = acc.ux.uxuxuy - acc.rho.uxuxuy * acc.in.ux;
            S.coeff(0, 10) = acc.ux.uxuyuy - acc.rho.uxuyuy * acc.in.ux;
            S.coeff(0, 11) = acc.ux.uyuyuy - acc.rho.uyuyuy * acc.in.ux;

            if constexpr (Rec)
            {
                S.coeff(0, 12) = acc.ux.uxmxx - acc.rho.uxmxx * acc.in.ux;
                S.coeff(0, 13) = acc.ux.uymxx - acc.rho.uymxx * acc.in.ux;
                S.coeff(0, 14) = acc.ux.uxmxy - acc.rho.uxmxy * acc.in.ux;
                S.coeff(0, 15) = acc.ux.uymxy - acc.rho.uymxy * acc.in.ux;
                S.coeff(0, 16) = acc.ux.uxmyy - acc.rho.uxmyy * acc.in.ux;
                S.coeff(0, 17) = acc.ux.uymyy - acc.rho.uymyy * acc.in.ux;
            }
        }

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

        if constexpr (RegOrder >= 3)
        {
            S.coeff(1, 8) = acc.uy.uxuxux - acc.rho.uxuxux * acc.in.uy;
            S.coeff(1, 9) = acc.uy.uxuxuy - acc.rho.uxuxuy * acc.in.uy;
            S.coeff(1, 10) = acc.uy.uxuyuy - acc.rho.uxuyuy * acc.in.uy;
            S.coeff(1, 11) = acc.uy.uyuyuy - acc.rho.uyuyuy * acc.in.uy;

            if constexpr (Rec)
            {
                S.coeff(1, 12) = acc.uy.uxmxx - acc.rho.uxmxx * acc.in.uy;
                S.coeff(1, 13) = acc.uy.uymxx - acc.rho.uymxx * acc.in.uy;
                S.coeff(1, 14) = acc.uy.uxmxy - acc.rho.uxmxy * acc.in.uy;
                S.coeff(1, 15) = acc.uy.uymxy - acc.rho.uymxy * acc.in.uy;
                S.coeff(1, 16) = acc.uy.uxmyy - acc.rho.uxmyy * acc.in.uy;
                S.coeff(1, 17) = acc.uy.uymyy - acc.rho.uymyy * acc.in.uy;
            }
        }

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

        if constexpr (RegOrder >= 3)
        {
            S.coeff(2, 8) = acc.mxx.uxuxux - acc.rho.uxuxux * acc.in.mxx;
            S.coeff(2, 9) = acc.mxx.uxuxuy - acc.rho.uxuxuy * acc.in.mxx;
            S.coeff(2, 10) = acc.mxx.uxuyuy - acc.rho.uxuyuy * acc.in.mxx;
            S.coeff(2, 11) = acc.mxx.uyuyuy - acc.rho.uyuyuy * acc.in.mxx;

            if constexpr (Rec)
            {
                S.coeff(2, 12) = acc.mxx.uxmxx - acc.rho.uxmxx * acc.in.mxx;
                S.coeff(2, 13) = acc.mxx.uymxx - acc.rho.uymxx * acc.in.mxx;
                S.coeff(2, 14) = acc.mxx.uxmxy - acc.rho.uxmxy * acc.in.mxx;
                S.coeff(2, 15) = acc.mxx.uymxy - acc.rho.uymxy * acc.in.mxx;
                S.coeff(2, 16) = acc.mxx.uxmyy - acc.rho.uxmyy * acc.in.mxx;
                S.coeff(2, 17) = acc.mxx.uymyy - acc.rho.uymyy * acc.in.mxx;
            }
        }

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

        if constexpr (RegOrder >= 3)
        {
            S.coeff(3, 8) = acc.mxy.uxuxux - acc.rho.uxuxux * acc.in.mxy;
            S.coeff(3, 9) = acc.mxy.uxuxuy - acc.rho.uxuxuy * acc.in.mxy;
            S.coeff(3, 10) = acc.mxy.uxuyuy - acc.rho.uxuyuy * acc.in.mxy;
            S.coeff(3, 11) = acc.mxy.uyuyuy - acc.rho.uyuyuy * acc.in.mxy;

            if constexpr (Rec)
            {
                S.coeff(3, 12) = acc.mxy.uxmxx - acc.rho.uxmxx * acc.in.mxy;
                S.coeff(3, 13) = acc.mxy.uymxx - acc.rho.uymxx * acc.in.mxy;
                S.coeff(3, 14) = acc.mxy.uxmxy - acc.rho.uxmxy * acc.in.mxy;
                S.coeff(3, 15) = acc.mxy.uymxy - acc.rho.uymxy * acc.in.mxy;
                S.coeff(3, 16) = acc.mxy.uxmyy - acc.rho.uxmyy * acc.in.mxy;
                S.coeff(3, 17) = acc.mxy.uymyy - acc.rho.uymyy * acc.in.mxy;
            }
        }

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

        if constexpr (RegOrder >= 3)
        {
            S.coeff(4, 8) = acc.myy.uxuxux - acc.rho.uxuxux * acc.in.myy;
            S.coeff(4, 9) = acc.myy.uxuxuy - acc.rho.uxuxuy * acc.in.myy;
            S.coeff(4, 10) = acc.myy.uxuyuy - acc.rho.uxuyuy * acc.in.myy;
            S.coeff(4, 11) = acc.myy.uyuyuy - acc.rho.uyuyuy * acc.in.myy;

            if constexpr (Rec)
            {
                S.coeff(4, 12) = acc.myy.uxmxx - acc.rho.uxmxx * acc.in.myy;
                S.coeff(4, 13) = acc.myy.uymxx - acc.rho.uymxx * acc.in.myy;
                S.coeff(4, 14) = acc.myy.uxmxy - acc.rho.uxmxy * acc.in.myy;
                S.coeff(4, 15) = acc.myy.uymxy - acc.rho.uymxy * acc.in.myy;
                S.coeff(4, 16) = acc.myy.uxmyy - acc.rho.uxmyy * acc.in.myy;
                S.coeff(4, 17) = acc.myy.uymyy - acc.rho.uymyy * acc.in.myy;
            }
        }

        S.b[4] = acc.rho.rho * acc.in.myy - acc.myy.rho;
    }
}