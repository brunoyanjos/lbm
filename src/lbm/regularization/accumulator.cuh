#pragma once

#include <cstdint>
#include "../state/lbm_state.cuh"
#include "../moment/moment_basis_value.cuh"

struct PopReconstructionAccumulator
{
    const LBMState &S;
    int c;
    size_t idx;
    int icx;
    int icy;
    real_t &f;

    template <auto Id>
    __device__ __forceinline__ void operator()()
    {
        const real_t a = device_field<Id>(S)[c][idx];
        f += a * moment_basis_value<Id>(icx, icy);
    }
};