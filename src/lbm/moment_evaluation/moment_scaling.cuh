#pragma once

#include "../../core/types.cuh"
#include "../stencil_active.cuh"

__device__ __forceinline__ void scale_to_stored_basis(real_t &ux, real_t &uy,
                                                      real_t &mxx, real_t &mxy, real_t &myy)
{
    ux *= Stencil::as2;
    uy *= Stencil::as2;
    mxx *= real_t(0.5) * Stencil::as4;
    mxy *= Stencil::as4;
    myy *= real_t(0.5) * Stencil::as4;
}
