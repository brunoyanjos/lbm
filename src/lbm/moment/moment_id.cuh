#pragma once

enum class MomentId : int
{
    rho = 0,

    ux,
    uy,

    mxx,
    mxy,
    myy,

    mxxx,
    mxxy,
    mxyy,
    myyy,

    mxxyy,
};

enum class NonlinearMomentId : int
{
    uxux = 0,
    uxuy,
    uyuy,

    uxuxux,
    uxuxuy,
    uxuyuy,
    uyuyuy,

    uxmxx,
    uymxx,
    uxmxy,
    uymxy,
    uxmyy,
    uymyy
};
