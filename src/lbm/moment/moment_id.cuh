#pragma once

enum class MomentId : int
{
    rho = 0,

    ux,
    uy,

    mxx,
    mxy,
    myy
};

enum class NonlinearMomentId : int
{
    uxux = 0,
    uxuy,
    uyuy
};