#pragma once

#include "core/types.cuh"

#include "lbm/hermite/hermite.cuh"
#include "lbm/moment/moment_id.cuh"

#include "lbm/boundary/fluid/accumulator.cuh"
#include "lbm/boundary/common/factor.cuh"

namespace boundary::fluid
{
    template <int RegOrder, bool Rec>
    __device__ __forceinline__ void evaluate_incoming(Accumulator<RegOrder, Rec> &acc,
                                                      const real_t *__restrict__ pop, int i)
    {
        acc.in.rho += pop[i];
        acc.in.ux += pop[i] * hermite<MomentId::ux>(i);
        acc.in.uy += pop[i] * hermite<MomentId::uy>(i);
        acc.in.mxx += pop[i] * hermite<MomentId::mxx>(i);
        acc.in.mxy += pop[i] * hermite<MomentId::mxy>(i);
        acc.in.myy += pop[i] * hermite<MomentId::myy>(i);

        acc.ux.rho += moment_factor<MomentId::rho, MomentId::ux>(i);
        acc.ux.ux += moment_factor<MomentId::ux, MomentId::ux>(i);
        acc.ux.uy += moment_factor<MomentId::uy, MomentId::ux>(i);
        acc.ux.mxx += moment_factor<MomentId::mxx, MomentId::ux>(i);
        acc.ux.mxy += moment_factor<MomentId::mxy, MomentId::ux>(i);
        acc.ux.myy += moment_factor<MomentId::myy, MomentId::ux>(i);

        if constexpr (RegOrder >= 3)
        {
            acc.ux.uxuxux += moment_factor<MomentId::mxxx, MomentId::ux>(i);
            acc.ux.uxuxuy += moment_factor<MomentId::mxxy, MomentId::ux>(i);
            acc.ux.uxuyuy += moment_factor<MomentId::mxyy, MomentId::ux>(i);
            acc.ux.uyuyuy += moment_factor<MomentId::myyy, MomentId::ux>(i);
        }

        acc.uy.rho += moment_factor<MomentId::rho, MomentId::uy>(i);
        acc.uy.ux += moment_factor<MomentId::ux, MomentId::uy>(i);
        acc.uy.uy += moment_factor<MomentId::uy, MomentId::uy>(i);
        acc.uy.mxx += moment_factor<MomentId::mxx, MomentId::uy>(i);
        acc.uy.mxy += moment_factor<MomentId::mxy, MomentId::uy>(i);
        acc.uy.myy += moment_factor<MomentId::myy, MomentId::uy>(i);

        if constexpr (RegOrder >= 3)
        {
            acc.uy.uxuxux += moment_factor<MomentId::mxxx, MomentId::uy>(i);
            acc.uy.uxuxuy += moment_factor<MomentId::mxxy, MomentId::uy>(i);
            acc.uy.uxuyuy += moment_factor<MomentId::mxyy, MomentId::uy>(i);
            acc.uy.uyuyuy += moment_factor<MomentId::myyy, MomentId::uy>(i);
        }

        acc.mxx.rho += moment_factor<MomentId::rho, MomentId::mxx>(i);
        acc.mxx.ux += moment_factor<MomentId::ux, MomentId::mxx>(i);
        acc.mxx.uy += moment_factor<MomentId::uy, MomentId::mxx>(i);
        acc.mxx.mxx += moment_factor<MomentId::mxx, MomentId::mxx>(i);
        acc.mxx.mxy += moment_factor<MomentId::mxy, MomentId::mxx>(i);
        acc.mxx.myy += moment_factor<MomentId::myy, MomentId::mxx>(i);

        if constexpr (RegOrder >= 3)
        {
            acc.mxx.uxuxux += moment_factor<MomentId::mxxx, MomentId::mxx>(i);
            acc.mxx.uxuxuy += moment_factor<MomentId::mxxy, MomentId::mxx>(i);
            acc.mxx.uxuyuy += moment_factor<MomentId::mxyy, MomentId::mxx>(i);
            acc.mxx.uyuyuy += moment_factor<MomentId::myyy, MomentId::mxx>(i);
        }

        acc.mxy.rho += moment_factor<MomentId::rho, MomentId::mxy>(i);
        acc.mxy.ux += moment_factor<MomentId::ux, MomentId::mxy>(i);
        acc.mxy.uy += moment_factor<MomentId::uy, MomentId::mxy>(i);
        acc.mxy.mxx += moment_factor<MomentId::mxx, MomentId::mxy>(i);
        acc.mxy.mxy += moment_factor<MomentId::mxy, MomentId::mxy>(i);
        acc.mxy.myy += moment_factor<MomentId::myy, MomentId::mxy>(i);

        if constexpr (RegOrder >= 3)
        {
            acc.mxy.uxuxux += moment_factor<MomentId::mxxx, MomentId::mxy>(i);
            acc.mxy.uxuxuy += moment_factor<MomentId::mxxy, MomentId::mxy>(i);
            acc.mxy.uxuyuy += moment_factor<MomentId::mxyy, MomentId::mxy>(i);
            acc.mxy.uyuyuy += moment_factor<MomentId::myyy, MomentId::mxy>(i);
        }

        acc.myy.rho += moment_factor<MomentId::rho, MomentId::myy>(i);
        acc.myy.ux += moment_factor<MomentId::ux, MomentId::myy>(i);
        acc.myy.uy += moment_factor<MomentId::uy, MomentId::myy>(i);
        acc.myy.mxx += moment_factor<MomentId::mxx, MomentId::myy>(i);
        acc.myy.mxy += moment_factor<MomentId::mxy, MomentId::myy>(i);
        acc.myy.myy += moment_factor<MomentId::myy, MomentId::myy>(i);

        if constexpr (RegOrder >= 3)
        {
            acc.myy.uxuxux += moment_factor<MomentId::mxxx, MomentId::myy>(i);
            acc.myy.uxuxuy += moment_factor<MomentId::mxxy, MomentId::myy>(i);
            acc.myy.uxuyuy += moment_factor<MomentId::mxyy, MomentId::myy>(i);
            acc.myy.uyuyuy += moment_factor<MomentId::myyy, MomentId::myy>(i);
        }
    }

    template <>
    __device__ __forceinline__ void evaluate_incoming(Accumulator<3, true> &acc,
                                                      const real_t *__restrict__ pop, int i)
    {
        acc.in.rho += pop[i];
        acc.in.ux += pop[i] * hermite<MomentId::ux>(i);
        acc.in.uy += pop[i] * hermite<MomentId::uy>(i);
        acc.in.mxx += pop[i] * hermite<MomentId::mxx>(i);
        acc.in.mxy += pop[i] * hermite<MomentId::mxy>(i);
        acc.in.myy += pop[i] * hermite<MomentId::myy>(i);

        acc.ux.rho += moment_factor<MomentId::rho, MomentId::ux>(i);
        acc.ux.ux += moment_factor<MomentId::ux, MomentId::ux>(i);
        acc.ux.uy += moment_factor<MomentId::uy, MomentId::ux>(i);
        acc.ux.mxx += moment_factor<MomentId::mxx, MomentId::ux>(i);
        acc.ux.mxy += moment_factor<MomentId::mxy, MomentId::ux>(i);
        acc.ux.myy += moment_factor<MomentId::myy, MomentId::ux>(i);
        acc.ux.uxuxux -= r::two * moment_factor<MomentId::mxxx, MomentId::ux>(i);
        acc.ux.uxuxuy -= r::two * moment_factor<MomentId::mxxy, MomentId::ux>(i);
        acc.ux.uxuyuy -= r::two * moment_factor<MomentId::mxyy, MomentId::ux>(i);
        acc.ux.uyuyuy -= r::two * moment_factor<MomentId::myyy, MomentId::ux>(i);
        acc.ux.uxmxx += r::three * moment_factor<MomentId::mxxx, MomentId::ux>(i);
        acc.ux.uymxx += moment_factor<MomentId::mxxy, MomentId::ux>(i);
        acc.ux.uxmxy += r::two * moment_factor<MomentId::mxxy, MomentId::ux>(i);
        acc.ux.uymxy += r::two * moment_factor<MomentId::mxyy, MomentId::ux>(i);
        acc.ux.uxmyy += moment_factor<MomentId::mxyy, MomentId::ux>(i);
        acc.ux.uymyy += r::three * moment_factor<MomentId::myyy, MomentId::ux>(i);

        acc.uy.rho += moment_factor<MomentId::rho, MomentId::uy>(i);
        acc.uy.ux += moment_factor<MomentId::ux, MomentId::uy>(i);
        acc.uy.uy += moment_factor<MomentId::uy, MomentId::uy>(i);
        acc.uy.mxx += moment_factor<MomentId::mxx, MomentId::uy>(i);
        acc.uy.mxy += moment_factor<MomentId::mxy, MomentId::uy>(i);
        acc.uy.myy += moment_factor<MomentId::myy, MomentId::uy>(i);
        acc.uy.uxuxux -= r::two * moment_factor<MomentId::mxxx, MomentId::uy>(i);
        acc.uy.uxuxuy -= r::two * moment_factor<MomentId::mxxy, MomentId::uy>(i);
        acc.uy.uxuyuy -= r::two * moment_factor<MomentId::mxyy, MomentId::uy>(i);
        acc.uy.uyuyuy -= r::two * moment_factor<MomentId::myyy, MomentId::uy>(i);
        acc.uy.uxmxx += r::three * moment_factor<MomentId::mxxx, MomentId::uy>(i);
        acc.uy.uymxx += moment_factor<MomentId::mxxy, MomentId::uy>(i);
        acc.uy.uxmxy += r::two * moment_factor<MomentId::mxxy, MomentId::uy>(i);
        acc.uy.uymxy += r::two * moment_factor<MomentId::mxyy, MomentId::uy>(i);
        acc.uy.uxmyy += moment_factor<MomentId::mxyy, MomentId::uy>(i);
        acc.uy.uymyy += r::three * moment_factor<MomentId::myyy, MomentId::uy>(i);

        acc.mxx.rho += moment_factor<MomentId::rho, MomentId::mxx>(i);
        acc.mxx.ux += moment_factor<MomentId::ux, MomentId::mxx>(i);
        acc.mxx.uy += moment_factor<MomentId::uy, MomentId::mxx>(i);
        acc.mxx.mxx += moment_factor<MomentId::mxx, MomentId::mxx>(i);
        acc.mxx.mxy += moment_factor<MomentId::mxy, MomentId::mxx>(i);
        acc.mxx.myy += moment_factor<MomentId::myy, MomentId::mxx>(i);
        acc.mxx.uxuxux -= r::two * moment_factor<MomentId::mxxx, MomentId::mxx>(i);
        acc.mxx.uxuxuy -= r::two * moment_factor<MomentId::mxxy, MomentId::mxx>(i);
        acc.mxx.uxuyuy -= r::two * moment_factor<MomentId::mxyy, MomentId::mxx>(i);
        acc.mxx.uyuyuy -= r::two * moment_factor<MomentId::myyy, MomentId::mxx>(i);
        acc.mxx.uxmxx += r::three * moment_factor<MomentId::mxxx, MomentId::mxx>(i);
        acc.mxx.uymxx += moment_factor<MomentId::mxxy, MomentId::mxx>(i);
        acc.mxx.uxmxy += r::two * moment_factor<MomentId::mxxy, MomentId::mxx>(i);
        acc.mxx.uymxy += r::two * moment_factor<MomentId::mxyy, MomentId::mxx>(i);
        acc.mxx.uxmyy += moment_factor<MomentId::mxyy, MomentId::mxx>(i);
        acc.mxx.uymyy += r::three * moment_factor<MomentId::myyy, MomentId::mxx>(i);

        acc.mxy.rho += moment_factor<MomentId::rho, MomentId::mxy>(i);
        acc.mxy.ux += moment_factor<MomentId::ux, MomentId::mxy>(i);
        acc.mxy.uy += moment_factor<MomentId::uy, MomentId::mxy>(i);
        acc.mxy.mxx += moment_factor<MomentId::mxx, MomentId::mxy>(i);
        acc.mxy.mxy += moment_factor<MomentId::mxy, MomentId::mxy>(i);
        acc.mxy.myy += moment_factor<MomentId::myy, MomentId::mxy>(i);
        acc.mxy.uxuxux -= r::two * moment_factor<MomentId::mxxx, MomentId::mxy>(i);
        acc.mxy.uxuxuy -= r::two * moment_factor<MomentId::mxxy, MomentId::mxy>(i);
        acc.mxy.uxuyuy -= r::two * moment_factor<MomentId::mxyy, MomentId::mxy>(i);
        acc.mxy.uyuyuy -= r::two * moment_factor<MomentId::myyy, MomentId::mxy>(i);
        acc.mxy.uxmxx += r::three * moment_factor<MomentId::mxxx, MomentId::mxy>(i);
        acc.mxy.uymxx += moment_factor<MomentId::mxxy, MomentId::mxy>(i);
        acc.mxy.uxmxy += r::two * moment_factor<MomentId::mxxy, MomentId::mxy>(i);
        acc.mxy.uymxy += r::two * moment_factor<MomentId::mxyy, MomentId::mxy>(i);
        acc.mxy.uxmyy += moment_factor<MomentId::mxyy, MomentId::mxy>(i);
        acc.mxy.uymyy += r::three * moment_factor<MomentId::myyy, MomentId::mxy>(i);

        acc.myy.rho += moment_factor<MomentId::rho, MomentId::myy>(i);
        acc.myy.ux += moment_factor<MomentId::ux, MomentId::myy>(i);
        acc.myy.uy += moment_factor<MomentId::uy, MomentId::myy>(i);
        acc.myy.mxx += moment_factor<MomentId::mxx, MomentId::myy>(i);
        acc.myy.mxy += moment_factor<MomentId::mxy, MomentId::myy>(i);
        acc.myy.myy += moment_factor<MomentId::myy, MomentId::myy>(i);
        acc.myy.uxuxux -= r::two * moment_factor<MomentId::mxxx, MomentId::myy>(i);
        acc.myy.uxuxuy -= r::two * moment_factor<MomentId::mxxy, MomentId::myy>(i);
        acc.myy.uxuyuy -= r::two * moment_factor<MomentId::mxyy, MomentId::myy>(i);
        acc.myy.uyuyuy -= r::two * moment_factor<MomentId::myyy, MomentId::myy>(i);
        acc.myy.uxmxx += r::three * moment_factor<MomentId::mxxx, MomentId::myy>(i);
        acc.myy.uymxx += moment_factor<MomentId::mxxy, MomentId::myy>(i);
        acc.myy.uxmxy += r::two * moment_factor<MomentId::mxxy, MomentId::myy>(i);
        acc.myy.uymxy += r::two * moment_factor<MomentId::mxyy, MomentId::myy>(i);
        acc.myy.uxmyy += moment_factor<MomentId::mxyy, MomentId::myy>(i);
        acc.myy.uymyy += r::three * moment_factor<MomentId::myyy, MomentId::myy>(i);
    }
}