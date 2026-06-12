#pragma once

#include <cstdint>
#include <cstdio>

#include "core/types.cuh"
#include "core/simulation_config.h"

#include "lbm/domain/mask_utils.cuh"
#include "lbm/stencil_active.cuh"
#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/system_data.cuh"
#include "lbm/boundary/common/gauss_elimination.cuh"

#include "lbm/boundary/dirichlet/accumulator.cuh"
#include "lbm/boundary/dirichlet/eval_density.cuh"
#include "lbm/boundary/dirichlet/evaluate_incoming.cuh"
#include "lbm/boundary/dirichlet/evaluate_outgoing.cuh"

namespace boundary::dirichlet
{

    template <int RegOrder, bool Rec, bool HighOrder>
    __device__ __forceinline__ void apply_boundary_numeric(
        real_t *__restrict__ pop,
        mask_t valid_mask,
        NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
    {
        const mask_t outgoing_mask = valid_mask;
        const mask_t incoming_mask = mask_opp(valid_mask);

        Accumulator acc{};

#pragma unroll
        for (int i = 0; i < Stencil::Q; ++i)
        {
            if (dir_valid(incoming_mask, i))
                evaluate_incoming(acc, M, pop, i);

            if (dir_valid(outgoing_mask, i))
                evaluate_outgoing(acc, M, i);
        }

        acc.in.normalize();

        M.mxy = (acc.rho.constant * acc.in.mxy - acc.mxy.constant) / (acc.mxy.mxy - acc.rho.mxy * acc.in.mxy);

        M.mxx = M.ux * M.ux;
        M.myy = M.uy * M.uy;

        M.rho = eval_density(acc, M);
    }

    template <bool HighOrder>
    __device__ __forceinline__ void apply_boundary_numeric(
        real_t *__restrict__ pop,
        mask_t valid_mask,
        NodeMomentsFor<3, false, HighOrder> &M)
    {
        const mask_t outgoing_mask = valid_mask;
        const mask_t incoming_mask = mask_opp(valid_mask);

        Accumulator acc{};

#pragma unroll
        for (int i = 0; i < Stencil::Q; ++i)
        {
            if (dir_valid(incoming_mask, i))
                evaluate_incoming(acc, M, pop, i);

            if (dir_valid(outgoing_mask, i))
                evaluate_outgoing(acc, M, i);
        }

        acc.in.normalize();

        M.mxy = (acc.rho.constant * acc.in.mxy - acc.mxy.constant) / (acc.mxy.mxy - acc.rho.mxy * acc.in.mxy);

        M.mxx = M.ux * M.ux;
        M.myy = M.uy * M.uy;
        M.mxxy = M.ux * M.ux * M.uy;
        M.mxyy = M.ux * M.uy * M.uy;

        if constexpr (HighOrder)
        {
            M.mxxx = M.ux * M.ux * M.ux;
            M.myyy = M.uy * M.uy * M.uy;
        }

        M.rho = eval_density(acc, M);
    }

    template <int RegOrder, bool Rec, bool HighOrder>
    __device__ __forceinline__ void apply_boundary_symbolic(
        real_t *__restrict__ pop,
        uint8_t node_id,
        NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
    {
        switch (static_cast<NodeId>(node_id))
        {
        case NodeId::NORTH_EAST:
        {
            // north-east corner formula
            const real_t rho_I = pop[0] + pop[1] + pop[2] + pop[5];
            const real_t inv_rho_I = r::one / rho_I;

            const real_t mxy_I = pop[5] * inv_rho_I;

            M.mxx = M.ux * M.ux;
            M.myy = M.uy * M.uy;

            M.rho = static_cast<real_t>(36) * (rho_I - mxy_I * rho_I + mxy_I * rho_I * OMEGA) /
                    (static_cast<real_t>(24) - static_cast<real_t>(18) * M.ux -
                     static_cast<real_t>(18) * M.ux * M.ux + OMEGA + static_cast<real_t>(3) * M.ux * OMEGA +
                     static_cast<real_t>(3) * M.ux * M.ux * OMEGA);

            M.mxy = (static_cast<real_t>(36) * mxy_I * rho_I - M.rho -
                     static_cast<real_t>(3) * M.ux * M.rho - static_cast<real_t>(3) * M.ux * M.ux * M.rho) /
                    (static_cast<real_t>(9) * M.rho);

            break;
        }
        case NodeId::NORTH_WEST:
        {
            // north-west corner formula
            const real_t rho_I = pop[0] + pop[2] + pop[3] + pop[6];
            const real_t inv_rho_I = r::one / rho_I;

            const real_t mxy_I = -pop[6] * inv_rho_I;

            M.mxx = M.ux * M.ux;
            M.myy = M.uy * M.uy;

            M.rho = -static_cast<real_t>(36) * (-rho_I - mxy_I * rho_I + mxy_I * rho_I * OMEGA) /
                    (static_cast<real_t>(24) + static_cast<real_t>(18) * M.ux -
                     static_cast<real_t>(18) * M.ux * M.ux + OMEGA - static_cast<real_t>(3) * M.ux * OMEGA +
                     static_cast<real_t>(3) * M.ux * M.ux * OMEGA);

            M.mxy = (static_cast<real_t>(36) * mxy_I * rho_I + M.rho -
                     static_cast<real_t>(3) * M.ux * M.rho + static_cast<real_t>(3) * M.ux * M.ux * M.rho) /
                    (static_cast<real_t>(9) * M.rho);

            break;
        }
        case NodeId::NORTH:
        {
            // north symbolic formula
            const real_t rho_I = pop[0] + pop[1] + pop[2] + pop[3] + pop[5] + pop[6];
            const real_t inv_rho_I = r::one / rho_I;

            const real_t mxy_I = (pop[5] - pop[6]) * inv_rho_I;

            M.mxx = M.ux * M.ux;
            M.myy = M.uy * M.uy;

            M.rho = static_cast<real_t>(6) * rho_I / static_cast<real_t>(5);
            M.mxy = (static_cast<real_t>(6) * mxy_I * rho_I - M.ux * M.rho) / (static_cast<real_t>(3) * M.rho);

            break;
        }
        case NodeId::EAST:
        {
            // east symbolic formula
            const real_t rho_I = pop[0] + pop[1] + pop[2] + pop[4] + pop[5] + pop[8];
            const real_t inv_rho_I = r::one / rho_I;

            const real_t mxy_I = (pop[5] - pop[8]) * inv_rho_I;

            M.mxx = M.ux * M.ux;
            M.myy = M.uy * M.uy;

            M.rho = static_cast<real_t>(6) * rho_I / static_cast<real_t>(5);
            M.mxy = static_cast<real_t>(2) * mxy_I * rho_I / M.rho;

            break;
        }
        case NodeId::WEST:
        {
            // west symbolic formula
            const real_t rho_I = pop[0] + pop[2] + pop[3] + pop[4] + pop[6] + pop[7];
            const real_t inv_rho_I = r::one / rho_I;

            const real_t mxy_I = (-pop[6] + pop[7]) * inv_rho_I;

            M.mxx = M.ux * M.ux;
            M.myy = M.uy * M.uy;

            M.rho = static_cast<real_t>(6) * rho_I / static_cast<real_t>(5);
            M.mxy = static_cast<real_t>(2) * mxy_I * rho_I / M.rho;

            break;
        }
        case NodeId::SOUTH:
        {
            // south symbolic formula
            const real_t rho_I = pop[0] + pop[1] + pop[3] + pop[4] + pop[7] + pop[8];
            const real_t inv_rho_I = r::one / rho_I;

            const real_t mxy_I = (pop[7] - pop[8]) * inv_rho_I;

            M.mxx = M.ux * M.ux;
            M.myy = M.uy * M.uy;

            M.rho = static_cast<real_t>(6) * rho_I / static_cast<real_t>(5);
            M.mxy = static_cast<real_t>(2) * mxy_I * rho_I / M.rho;

            break;
        }
        case NodeId::SOUTH_EAST:
        {
            // south-east corner formula
            const real_t rho_I = pop[0] + pop[1] + pop[4] + pop[8];
            const real_t inv_rho_I = r::one / rho_I;

            const real_t mxy_I = -pop[8] * inv_rho_I;

            M.mxx = M.ux * M.ux;
            M.mxy = M.ux * M.uy;
            M.myy = M.uy * M.uy;

            M.rho = -static_cast<real_t>(36) * (-rho_I - mxy_I * rho_I + mxy_I * rho_I * OMEGA) /
                    (static_cast<real_t>(24) + OMEGA);
            M.mxy = (static_cast<real_t>(36) * mxy_I * rho_I + M.rho) / (static_cast<real_t>(24) + OMEGA);

            break;
        }
        case NodeId::SOUTH_WEST:
        {
            // south-west corner formula
            const real_t rho_I = pop[0] + pop[3] + pop[4] + pop[7];
            const real_t inv_rho_I = r::one / rho_I;

            const real_t mxy_I = pop[7] * inv_rho_I;

            M.mxx = M.ux * M.ux;
            M.myy = M.uy * M.uy;

            M.rho = static_cast<real_t>(36) * (rho_I - mxy_I * rho_I + mxy_I * rho_I * OMEGA) /
                    (static_cast<real_t>(24) + OMEGA);
            M.mxy = (static_cast<real_t>(36) * mxy_I * rho_I - M.rho) / (static_cast<real_t>(24) + OMEGA);

            break;
        }
        default:
            // optional fallback
            break;
        }
    }

    template <bool HighOrder>
    __device__ __forceinline__ void apply_boundary_symbolic(
        real_t *__restrict__ pop,
        uint8_t node_id,
        NodeMomentsFor<3, false, HighOrder> &M)
    {
        (void)pop;
        (void)node_id;
        (void)M;
    }

    template <int RegOrder, bool Rec, bool HighOrder>
    __device__ __forceinline__ void apply_boundary(
        real_t *__restrict__ pop,
        uint8_t node_id,
        mask_t valid_mask,
        NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
    {
        if constexpr (USE_SYMBOLIC_BOUNDARY)
        {
            apply_boundary_symbolic(pop, node_id, M);
        }
        else
        {
            apply_boundary_numeric(pop, valid_mask, M);
        }
    }

    template <bool HighOrder>
    __device__ __forceinline__ void apply_boundary(
        real_t *__restrict__ pop,
        uint8_t node_id,
        mask_t valid_mask,
        NodeMomentsFor<3, false, HighOrder> &M)
    {
        if constexpr (USE_SYMBOLIC_BOUNDARY)
        {
            apply_boundary_symbolic(pop, node_id, M);
        }
        else
        {
            apply_boundary_numeric(pop, valid_mask, M);
        }
    }

    template <int RegOrder, bool Rec, bool HighOrder>
    __device__ __forceinline__ void apply_boundary(
        real_t *__restrict__ pop,
        mask_t valid_mask,
        NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
    {
        apply_boundary_numeric(pop, valid_mask, M);
    }

    template <bool HighOrder>
    __device__ __forceinline__ void apply_boundary(
        real_t *__restrict__ pop,
        mask_t valid_mask,
        NodeMomentsFor<3, false, HighOrder> &M)
    {
        apply_boundary_numeric(pop, valid_mask, M);
    }

}
