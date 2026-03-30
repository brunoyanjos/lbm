#pragma once

#include <cstdint>
#include <cstdio>

#include "core/types.cuh"

#include "lbm/domain/mask_utils.cuh"
#include "lbm/domain/domain_tags.cuh"
#include "lbm/stencil_active.cuh"
#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/system_data.cuh"
#include "lbm/boundary/common/gauss_elimination.cuh"
#include "lbm/boundary/common/layout/unknown_moments.cuh"

#include "lbm/boundary/dirichlet/build_boundary_system.cuh"
#include "lbm/boundary/dirichlet/dirichlet_accumulator.cuh"
#include "lbm/boundary/dirichlet/dirichlet_eval_density.cuh"
#include "lbm/boundary/dirichlet/dirichlet_incoming_evaluation.cuh"
#include "lbm/boundary/dirichlet/dirichlet_unknown_moments.cuh"
#include "lbm/boundary/dirichlet/dirichlet_outgoing_evaluation.cuh"

__device__ __forceinline__ void apply_boundary(
    real_t *__restrict__ pop,
    uint8_t node_id,
    NodeMoments &M)
{
    switch (node_id)
    {
    case to_u8(NodeId::ONE):
    {
        const real_t rhoI = pop[0] + pop[1] + pop[2] + pop[5];
        const real_t inv_rhoI = r::one / rhoI;

        const real_t mxxI = (pop[1] + pop[5]) * inv_rhoI - Stencil::cs2;
        const real_t mxyI = (pop[5]) * inv_rhoI;
        const real_t myyI = (pop[2] + pop[5]) * inv_rhoI - Stencil::cs2;

        M.rho =
            (real_t(12) *
             (-real_t(3) * rhoI - real_t(3) * mxxI * rhoI + real_t(7) * mxyI * rhoI - real_t(3) * myyI * rhoI + real_t(3) * mxxI * rhoI * OMEGA - real_t(7) * mxyI * rhoI * OMEGA + real_t(3) * myyI * rhoI * OMEGA)) /
            (-real_t(16) + real_t(14) * M.ux + real_t(14) * M.uy - real_t(9) * OMEGA + M.ux * OMEGA + real_t(15) * M.ux * M.ux * OMEGA + M.uy * OMEGA - real_t(9) * M.ux * M.uy * OMEGA + real_t(15) * M.uy * M.uy * OMEGA);

        M.mxx =
            -(real_t(2) *
              (-real_t(9) * mxxI * rhoI + real_t(6) * mxyI * rhoI - M.rho + real_t(2) * M.ux * M.rho - M.uy * M.rho)) /
            (real_t(9) * M.rho);

        M.mxy =
            -(real_t(18) * mxxI * rhoI - real_t(132) * mxyI * rhoI + real_t(18) * myyI * rhoI + real_t(7) * M.rho + real_t(7) * M.ux * M.rho + real_t(7) * M.uy * M.rho) /
            (real_t(27) * M.rho);

        M.myy =
            (real_t(2) *
             (-real_t(6) * mxyI * rhoI + real_t(9) * myyI * rhoI + M.rho + M.ux * M.rho - real_t(2) * M.uy * M.rho)) /
            (real_t(9) * M.rho);

        break;
    }
    case to_u8(NodeId::TWO):
    {
        const real_t rhoI = pop[0] + pop[2] + pop[3] + pop[6];
        const real_t inv_rhoI = r::one / rhoI;

        const real_t mxxI = (pop[3] + pop[6]) * inv_rhoI - Stencil::cs2;
        const real_t mxyI = (-pop[6]) * inv_rhoI;
        const real_t myyI = (pop[2] + pop[6]) * inv_rhoI - Stencil::cs2;

        M.rho =
            (real_t(12) * rhoI *
             (-real_t(3) - real_t(3) * mxxI - real_t(7) * mxyI - real_t(3) * myyI + real_t(3) * mxxI * OMEGA + real_t(7) * mxyI * OMEGA + real_t(3) * myyI * OMEGA)) /
            (-real_t(16) - real_t(14) * M.ux + real_t(14) * M.uy - real_t(9) * OMEGA - M.ux * OMEGA + real_t(15) * M.ux * M.ux * OMEGA + M.uy * OMEGA + real_t(9) * M.ux * M.uy * OMEGA + real_t(15) * M.uy * M.uy * OMEGA);

        M.mxx =
            (real_t(2) *
             (real_t(9) * mxxI * rhoI + real_t(6) * mxyI * rhoI + M.rho + real_t(2) * M.ux * M.rho + M.uy * M.rho)) /
            (real_t(9) * M.rho);

        M.mxy =
            -(-real_t(18) * mxxI * rhoI - real_t(132) * mxyI * rhoI - real_t(18) * myyI * rhoI - real_t(7) * M.rho + real_t(7) * M.ux * M.rho - real_t(7) * M.uy * M.rho) /
            (real_t(27) * M.rho);

        M.myy =
            -(real_t(2) *
              (-real_t(6) * mxyI * rhoI - real_t(9) * myyI * rhoI - M.rho + M.ux * M.rho + real_t(2) * M.uy * M.rho)) /
            (real_t(9) * M.rho);

        break;
    }
    case to_u8(NodeId::THREE):
    {
        const real_t rhoI = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];
        const real_t inv_rhoI = r::one / rhoI;

        const real_t mxxI = (pop[1] + pop[3] + pop[5] + pop[6]) * inv_rhoI - Stencil::cs2;
        const real_t mxyI = (pop[5] - pop[6]) * inv_rhoI;
        const real_t myyI = (pop[2] + pop[5] + pop[6]) * inv_rhoI - Stencil::cs2;

        M.rho =
            (real_t(3) *
             (-real_t(4) * rhoI - real_t(3) * myyI * rhoI + real_t(3) * myyI * rhoI * OMEGA)) /
            (-real_t(9) + real_t(3) * M.uy - OMEGA + real_t(3) * M.uy * OMEGA + real_t(6) * M.uy * M.uy * OMEGA);

        M.mxx =
            (real_t(6) * mxxI * rhoI) /
            (real_t(5) * M.rho);

        M.mxy =
            -(-real_t(6) * mxyI * rhoI + M.ux * M.rho) /
            (real_t(3) * M.rho);

        M.myy =
            -(-real_t(9) * myyI * rhoI - M.rho + real_t(3) * M.uy * M.rho) /
            (real_t(6) * M.rho);

        break;
    }
    case to_u8(NodeId::FOUR):
    {
        const real_t rhoI = pop[0] + pop[1] + pop[4] + pop[8];
        const real_t inv_rhoI = r::one / rhoI;

        const real_t mxxI = (pop[1] + pop[8]) * inv_rhoI - Stencil::cs2;
        const real_t mxyI = (-pop[8]) * inv_rhoI;
        const real_t myyI = (pop[4] + pop[8]) * inv_rhoI - Stencil::cs2;

        M.rho =
            (real_t(12) * rhoI *
             (-real_t(3) - real_t(3) * mxxI - real_t(7) * mxyI - real_t(3) * myyI + real_t(3) * mxxI * OMEGA + real_t(7) * mxyI * OMEGA + real_t(3) * myyI * OMEGA)) /
            (-real_t(16) + real_t(14) * M.ux - real_t(14) * M.uy - real_t(9) * OMEGA + M.ux * OMEGA + real_t(15) * M.ux * M.ux * OMEGA - M.uy * OMEGA + real_t(9) * M.ux * M.uy * OMEGA + real_t(15) * M.uy * M.uy * OMEGA);

        M.mxx =
            -(real_t(2) *
              (-real_t(9) * mxxI * rhoI - real_t(6) * mxyI * rhoI - M.rho + real_t(2) * M.ux * M.rho + M.uy * M.rho)) /
            (real_t(9) * M.rho);

        M.mxy =
            -(-real_t(18) * mxxI * rhoI - real_t(132) * mxyI * rhoI - real_t(18) * myyI * rhoI - real_t(7) * M.rho - real_t(7) * M.ux * M.rho + real_t(7) * M.uy * M.rho) /
            (real_t(27) * M.rho);

        M.myy =
            (real_t(2) *
             (real_t(6) * mxyI * rhoI + real_t(9) * myyI * rhoI + M.rho + M.ux * M.rho + real_t(2) * M.uy * M.rho)) /
            (real_t(9) * M.rho);

        break;
    }
    case to_u8(NodeId::FIVE):
    {
        const real_t rhoI = pop[0] + pop[1] + pop[2] + pop[4] + pop[5] + pop[8];
        const real_t inv_rhoI = r::one / rhoI;

        const real_t mxxI = (pop[1] + pop[5] + pop[8]) * inv_rhoI - Stencil::cs2;
        const real_t mxyI = (pop[5] - pop[8]) * inv_rhoI;
        const real_t myyI = (pop[2] + pop[4] + pop[5] + pop[8]) * inv_rhoI - Stencil::cs2;

        M.rho =
            (real_t(3) *
             (-real_t(4) * rhoI - real_t(3) * mxxI * rhoI + real_t(3) * mxxI * rhoI * OMEGA)) /
            (-real_t(9) + real_t(3) * M.ux - OMEGA + real_t(3) * M.ux * OMEGA + real_t(6) * M.ux * M.ux * OMEGA);

        M.mxx =
            -(-real_t(9) * mxxI * rhoI - M.rho + real_t(3) * M.ux * M.rho) /
            (real_t(6) * M.rho);

        M.mxy =
            -(-real_t(6) * mxyI * rhoI + M.uy * M.rho) /
            (real_t(3) * M.rho);

        M.myy =
            (real_t(6) * myyI * rhoI) /
            (real_t(5) * M.rho);

        break;
    }
    case to_u8(NodeId::SEVEN):
    {
        const real_t rhoI = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[8];
        const real_t inv_rhoI = r::one / rhoI;

        const real_t mxxI = (pop[1] + pop[3] + pop[5] + pop[6] + pop[8]) * inv_rhoI - Stencil::cs2;
        const real_t mxyI = (pop[5] - pop[6] - pop[8]) * inv_rhoI;
        const real_t myyI = (pop[2] + pop[4] + pop[5] + pop[6] + pop[8]) * inv_rhoI - Stencil::cs2;

        M.rho =
            (real_t(36) * rhoI *
             (-real_t(23) - real_t(3) * mxxI - real_t(9) * mxyI - real_t(3) * myyI + real_t(3) * mxxI * OMEGA + real_t(9) * mxyI * OMEGA + real_t(3) * myyI * OMEGA)) /
            (-real_t(792) + real_t(30) * M.ux + real_t(30) * M.uy - real_t(13) * OMEGA + real_t(39) * M.ux * OMEGA + real_t(69) * M.ux * M.ux * OMEGA + real_t(39) * M.uy * OMEGA + real_t(207) * M.ux * M.uy * OMEGA + real_t(69) * M.uy * M.uy * OMEGA);

        M.mxx =
            -(-real_t(75) * mxxI * rhoI - real_t(18) * mxyI * rhoI - real_t(6) * myyI * rhoI - real_t(2) * M.rho + real_t(6) * M.ux * M.rho + real_t(6) * M.uy * M.rho) /
            (real_t(69) * M.rho);

        M.mxy =
            -(-real_t(3) * mxxI * rhoI - real_t(32) * mxyI * rhoI - real_t(3) * myyI * rhoI - M.rho + real_t(3) * M.ux * M.rho + real_t(3) * M.uy * M.rho) /
            (real_t(23) * M.rho);

        M.myy =
            -(-real_t(6) * mxxI * rhoI - real_t(18) * mxyI * rhoI - real_t(75) * myyI * rhoI - real_t(2) * M.rho + real_t(6) * M.ux * M.rho + real_t(6) * M.uy * M.rho) /
            (real_t(69) * M.rho);

        break;
    }
    case to_u8(NodeId::EIGHT):
    {
        const real_t rhoI = pop[0] + pop[3] + pop[4] + pop[7];
        const real_t inv_rhoI = r::one / rhoI;

        const real_t mxxI = (pop[3] + pop[7]) * inv_rhoI - Stencil::cs2;
        const real_t mxyI = (pop[7]) * inv_rhoI;
        const real_t myyI = (pop[4] + pop[7]) * inv_rhoI - Stencil::cs2;

        M.rho =
            (real_t(12) * rhoI *
             (-real_t(3) - real_t(3) * mxxI + real_t(7) * mxyI - real_t(3) * myyI + real_t(3) * mxxI * OMEGA - real_t(7) * mxyI * OMEGA + real_t(3) * myyI * OMEGA)) /
            (-real_t(16) - real_t(14) * M.ux - real_t(14) * M.uy - real_t(9) * OMEGA - M.ux * OMEGA + real_t(15) * M.ux * M.ux * OMEGA - M.uy * OMEGA - real_t(9) * M.ux * M.uy * OMEGA + real_t(15) * M.uy * M.uy * OMEGA);

        M.mxx =
            (real_t(2) *
             (real_t(9) * mxxI * rhoI - real_t(6) * mxyI * rhoI + M.rho + real_t(2) * M.ux * M.rho - M.uy * M.rho)) /
            (real_t(9) * M.rho);

        M.mxy =
            -(real_t(18) * mxxI * rhoI - real_t(132) * mxyI * rhoI + real_t(18) * myyI * rhoI + real_t(7) * M.rho - real_t(7) * M.ux * M.rho - real_t(7) * M.uy * M.rho) /
            (real_t(27) * M.rho);

        M.myy =
            -(real_t(2) *
              (real_t(6) * mxyI * rhoI - real_t(9) * myyI * rhoI - M.rho + M.ux * M.rho - real_t(2) * M.uy * M.rho)) /
            (real_t(9) * M.rho);

        break;
    }
    case to_u8(NodeId::TEN):
    {
        const real_t rhoI = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
        const real_t inv_rhoI = r::one / rhoI;

        const real_t mxxI = (pop[3] + pop[6] + pop[7]) * inv_rhoI - Stencil::cs2;
        const real_t mxyI = (-pop[6] + pop[7]) * inv_rhoI;
        const real_t myyI = (pop[2] + pop[4] + pop[6] + pop[7]) * inv_rhoI - Stencil::cs2;

        M.rho =
            (real_t(3) *
             (-real_t(4) * rhoI - real_t(3) * mxxI * rhoI + real_t(3) * mxxI * rhoI * OMEGA)) /
            (-real_t(9) - real_t(3) * M.ux - OMEGA - real_t(3) * M.ux * OMEGA + real_t(6) * M.ux * M.ux * OMEGA);

        M.mxx =
            -(-real_t(9) * mxxI * rhoI - M.rho - real_t(3) * M.ux * M.rho) /
            (real_t(6) * M.rho);

        M.mxy =
            -(-real_t(6) * mxyI * rhoI - M.uy * M.rho) /
            (real_t(3) * M.rho);

        M.myy =
            (real_t(6) * myyI * rhoI) /
            (real_t(5) * M.rho);

        break;
    }
    case to_u8(NodeId::ELEVEN):
    {
        const real_t rhoI = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[6] + pop[7];
        const real_t inv_rhoI = r::one / rhoI;

        const real_t mxxI = (pop[1] + pop[3] + pop[5] + pop[6] + pop[7]) * inv_rhoI - Stencil::cs2;
        const real_t mxyI = (pop[5] - pop[6] + pop[7]) * inv_rhoI;
        const real_t myyI = (pop[2] + pop[4] + pop[5] + pop[6] + pop[7]) * inv_rhoI - Stencil::cs2;

        M.rho =
            (real_t(36) * rhoI *
             (-real_t(23) - real_t(3) * mxxI + real_t(9) * mxyI - real_t(3) * myyI + real_t(3) * mxxI * OMEGA - real_t(9) * mxyI * OMEGA + real_t(3) * myyI * OMEGA)) /
            (-real_t(792) - real_t(30) * M.ux + real_t(30) * M.uy - real_t(13) * OMEGA - real_t(39) * M.ux * OMEGA + real_t(69) * M.ux * M.ux * OMEGA + real_t(39) * M.uy * OMEGA - real_t(207) * M.ux * M.uy * OMEGA + real_t(69) * M.uy * M.uy * OMEGA);

        M.mxx =
            -(-real_t(75) * mxxI * rhoI + real_t(18) * mxyI * rhoI - real_t(6) * myyI * rhoI - real_t(2) * M.rho - real_t(6) * M.ux * M.rho + real_t(6) * M.uy * M.rho) /
            (real_t(69) * M.rho);

        M.mxy =
            -(real_t(3) * mxxI * rhoI - real_t(32) * mxyI * rhoI + real_t(3) * myyI * rhoI + M.rho + real_t(3) * M.ux * M.rho - real_t(3) * M.uy * M.rho) /
            (real_t(23) * M.rho);

        M.myy =
            -(-real_t(6) * mxxI * rhoI + real_t(18) * mxyI * rhoI - real_t(75) * myyI * rhoI - real_t(2) * M.rho - real_t(6) * M.ux * M.rho + real_t(6) * M.uy * M.rho) /
            (real_t(69) * M.rho);

        break;
    }
    case to_u8(NodeId::TWELVE):
    {
        const real_t rhoI = pop[0] + pop[1] + pop[3] + pop[4] + pop[7] + pop[8];
        const real_t inv_rhoI = r::one / rhoI;

        const real_t mxxI = (pop[1] + pop[3] + pop[7] + pop[8]) * inv_rhoI - Stencil::cs2;
        const real_t mxyI = (pop[7] - pop[8]) * inv_rhoI;
        const real_t myyI = (pop[4] + pop[7] + pop[8]) * inv_rhoI - Stencil::cs2;

        M.rho =
            (real_t(3) *
             (-real_t(4) * rhoI - real_t(3) * myyI * rhoI + real_t(3) * myyI * rhoI * OMEGA)) /
            (-real_t(9) - real_t(3) * M.uy - OMEGA - real_t(3) * M.uy * OMEGA + real_t(6) * M.uy * M.uy * OMEGA);

        M.mxx =
            (real_t(6) * mxxI * rhoI) /
            (real_t(5) * M.rho);

        M.mxy =
            -(-real_t(6) * mxyI * rhoI - M.ux * M.rho) /
            (real_t(3) * M.rho);

        M.myy =
            -(-real_t(9) * myyI * rhoI - M.rho - real_t(3) * M.uy * M.rho) /
            (real_t(6) * M.rho);

        break;
    }
    case to_u8(NodeId::THIRTEEN):
    {
        const real_t rhoI = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[5] + pop[7] + pop[8];
        const real_t inv_rhoI = r::one / rhoI;

        const real_t mxxI = (pop[1] + pop[3] + pop[5] + pop[7] + pop[8]) * inv_rhoI - Stencil::cs2;
        const real_t mxyI = (pop[5] + pop[7] - pop[8]) * inv_rhoI;
        const real_t myyI = (pop[2] + pop[4] + pop[5] + pop[7] + pop[8]) * inv_rhoI - Stencil::cs2;

        M.rho =
            (real_t(36) * rhoI *
             (-real_t(23) - real_t(3) * mxxI + real_t(9) * mxyI - real_t(3) * myyI + real_t(3) * mxxI * OMEGA - real_t(9) * mxyI * OMEGA + real_t(3) * myyI * OMEGA)) /
            (-real_t(792) + real_t(30) * M.ux - real_t(30) * M.uy - real_t(13) * OMEGA + real_t(39) * M.ux * OMEGA + real_t(69) * M.ux * M.ux * OMEGA - real_t(39) * M.uy * OMEGA - real_t(207) * M.ux * M.uy * OMEGA + real_t(69) * M.uy * M.uy * OMEGA);

        M.mxx =
            -(-real_t(75) * mxxI * rhoI + real_t(18) * mxyI * rhoI - real_t(6) * myyI * rhoI - real_t(2) * M.rho + real_t(6) * M.ux * M.rho - real_t(6) * M.uy * M.rho) /
            (real_t(69) * M.rho);

        M.mxy =
            (real_t(3) * mxxI * rhoI - real_t(32) * mxyI * rhoI + real_t(3) * myyI * rhoI + M.rho - real_t(3) * M.ux * M.rho + real_t(3) * M.uy * M.rho) /
            (real_t(23) * M.rho);

        M.myy =
            -(-real_t(6) * mxxI * rhoI + real_t(18) * mxyI * rhoI - real_t(75) * myyI * rhoI - real_t(2) * M.rho + real_t(6) * M.ux * M.rho - real_t(6) * M.uy * M.rho) /
            (real_t(69) * M.rho);

        break;
    }
    case to_u8(NodeId::FOURTEEN):
    {
        const real_t rhoI = pop[0] + pop[1] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7] + pop[8];
        const real_t inv_rhoI = r::one / rhoI;

        const real_t mxxI = (pop[1] + pop[3] + pop[6] + pop[7] + pop[8]) * inv_rhoI - Stencil::cs2;
        const real_t mxyI = (-pop[6] + pop[7] - pop[8]) * inv_rhoI;
        const real_t myyI = (pop[2] + pop[4] + pop[6] + pop[7] + pop[8]) * inv_rhoI - Stencil::cs2;

        M.rho =
            (real_t(36) * rhoI *
             (-real_t(23) - real_t(3) * mxxI - real_t(9) * mxyI - real_t(3) * myyI + real_t(3) * mxxI * OMEGA + real_t(9) * mxyI * OMEGA + real_t(3) * myyI * OMEGA)) /
            (-real_t(792) - real_t(30) * M.ux - real_t(30) * M.uy - real_t(13) * OMEGA - real_t(39) * M.ux * OMEGA + real_t(69) * M.ux * M.ux * OMEGA - real_t(39) * M.uy * OMEGA + real_t(207) * M.ux * M.uy * OMEGA + real_t(69) * M.uy * M.uy * OMEGA);

        M.mxx =
            -(-real_t(75) * mxxI * rhoI - real_t(18) * mxyI * rhoI - real_t(6) * myyI * rhoI - real_t(2) * M.rho - real_t(6) * M.ux * M.rho - real_t(6) * M.uy * M.rho) /
            (real_t(69) * M.rho);

        M.mxy =
            -(-real_t(3) * mxxI * rhoI - real_t(32) * mxyI * rhoI - real_t(3) * myyI * rhoI - M.rho - real_t(3) * M.ux * M.rho - real_t(3) * M.uy * M.rho) /
            (real_t(23) * M.rho);

        M.myy =
            -(-real_t(6) * mxxI * rhoI - real_t(18) * mxyI * rhoI - real_t(75) * myyI * rhoI - real_t(2) * M.rho - real_t(6) * M.ux * M.rho - real_t(6) * M.uy * M.rho) /
            (real_t(69) * M.rho);

        break;
    }
    default:
    {
        break;
    }
    }
}