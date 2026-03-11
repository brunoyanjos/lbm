#include "lbm_mom_step.cuh"
#include "stencil_active.cuh"
#include "boundary/regularized_boundary_condition.cuh"
#include "boundary/fluid_moment_evaluation.cuh"
#include "boundary/fluid_boundary_evaluation.cuh"
#include "collision/collision.cuh"
#include "moment_evaluation/moment_evaluation.cuh"
#include "moment_evaluation/moment_scaling.cuh"
#include "state/state_store.cuh"
#include "population/pop_reconstruction.cuh"

#include "../geometries/active_geometry.cuh"

#include "../core/indexing.cuh"
#include "../core/cuda_utils.cuh"
#include "../core/cuda_utils.cuh"

__global__ void lbm_mom_step_kernel(LBMState S, DomainTags T)
{
    int x, y;
    const size_t idx = idxThreadGlobal2D(x, y);
    if (idx == INVALID_INDEX)
        return;

    const uint8_t wall_id = T.d_node[idx];
    const uint32_t valid_ms = T.d_valid[idx];

    if (wall_id == to_u8(NodeId::SOLID))
        return;

    real_t pop[Stencil::Q];

    // 1) reconstruct streamed populations from neighbor moments
    // reconstruct_streamed_pop(pop, S, c, x, y);

    real_t rho, ux, uy, mxx, mxy, myy;

    if (wall_id != static_cast<uint8_t>(NodeId::FLUID))
    {
        Geometry::bc_velocity(x, y, ux, uy);

        apply_boundary_dirichlet(pop, valid_ms, rho, ux, uy, mxx, mxy, myy);
    }
    else
    {
        if (is_full_mask(valid_ms))
        {
            evaluate_moments_from_pop(pop, rho, ux, uy, mxx, mxy, myy);
        }
        else
        {
            // rho = S.d_rho[c][idx] + Geometry::RHO_0;
            // ux = S.d_ux[c][idx] / Stencil::as2;
            // uy = S.d_uy[c][idx] / Stencil::as2;
            // mxx = S.d_mxx[c][idx] / (Stencil::as4 * r::half);
            // mxy = S.d_mxy[c][idx] / Stencil::as4;
            // myy = S.d_myy[c][idx] / (Stencil::as4 * r::half);

            evaluate_fluid_node(pop, valid_ms, rho, ux, uy, mxx, mxy, myy);
        }
    }

    // 3) scale to the stored basis
    scale_to_stored_basis(ux, uy, mxx, mxy, myy);

    // 4) collide in moment space
    moment_space_collision(ux, uy, mxx, mxy, myy);

    // 5) store next
    // store_next_state(S, n, idx, rho, ux, uy, mxx, mxy, myy);
}

__global__ void load_layers(LBMState S)
{
    const int x = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
                  static_cast<int>(threadIdx.x);
    if (x >= Geometry::NX)
        return;

    const size_t idx_0 = x + 0 * Geometry::NX;
    const size_t idx_1 = x + 1 * Geometry::NX;
    const size_t idx_2 = x + 2 * Geometry::NX;

    S.d_rho[idx_0] = S.global_rho[idx_0];
    S.d_ux[idx_0] = S.global_ux[idx_0];
    S.d_uy[idx_0] = S.global_uy[idx_0];
    S.d_mxx[idx_0] = S.global_mxx[idx_0];
    S.d_mxy[idx_0] = S.global_mxy[idx_0];
    S.d_myy[idx_0] = S.global_myy[idx_0];

    S.d_rho[idx_1] = S.global_rho[idx_1];
    S.d_ux[idx_1] = S.global_ux[idx_1];
    S.d_uy[idx_1] = S.global_uy[idx_1];
    S.d_mxx[idx_1] = S.global_mxx[idx_1];
    S.d_mxy[idx_1] = S.global_mxy[idx_1];
    S.d_myy[idx_1] = S.global_myy[idx_1];

    S.d_rho[idx_2] = S.global_rho[idx_2];
    S.d_ux[idx_2] = S.global_ux[idx_2];
    S.d_uy[idx_2] = S.global_uy[idx_2];
    S.d_mxx[idx_2] = S.global_mxx[idx_2];
    S.d_mxy[idx_2] = S.global_mxy[idx_2];
    S.d_myy[idx_2] = S.global_myy[idx_2];
}

__global__ void evaluate_first_layer(LBMState S, DomainTags T)
{
    const int x = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
                  static_cast<int>(threadIdx.x);

    if (x >= Geometry::NX)
        return;

    int idx = x + 0 * Geometry::NX;

    const uint8_t wall_id = T.d_node[idx];
    const uint32_t valid_ms = T.d_valid[idx];

    real_t pop[Stencil::Q];

    reconstruct_streamed_pop(pop, S, x, 0);

    real_t rho, ux, uy, mxx, mxy, myy;

    if (wall_id != to_u8(NodeId::FLUID))
    {
        Geometry::bc_velocity(x, 0, ux, uy);

        apply_boundary_dirichlet(pop, valid_ms, rho, ux, uy, mxx, mxy, myy);
    }

    scale_to_stored_basis(ux, uy, mxx, mxy, myy);
    moment_space_collision(ux, uy, mxx, mxy, myy);

    S.global_rho[idx] = rho - Geometry::RHO_0;
    S.global_ux[idx] = ux;
    S.global_uy[idx] = uy;
    S.global_mxx[idx] = mxx;
    S.global_mxy[idx] = mxy;
    S.global_myy[idx] = myy;
}

__global__ void evaluate_middle_layer(LBMState S, DomainTags T, int y0)
{
    const int x = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
                  static_cast<int>(threadIdx.x);

    if (x >= Geometry::NX)
        return;

    const size_t layer_idx = x + 1 * Geometry::NX;
    const size_t global_idx = x + (y0 + 1) * Geometry::NX;

    const uint8_t wall_id = T.d_node[global_idx];
    const uint32_t valid_ms = T.d_valid[global_idx];

    real_t pop[Stencil::Q];

    reconstruct_streamed_pop(pop, S, x, 1);

    real_t rho, ux, uy, mxx, mxy, myy;

    if (wall_id != to_u8(NodeId::FLUID))
    {
        Geometry::bc_velocity(x, y0 + 1, ux, uy);

        apply_boundary_dirichlet(pop, valid_ms, rho, ux, uy, mxx, mxy, myy);
    }
    else
    {
        evaluate_moments_from_pop(pop, rho, ux, uy, mxx, mxy, myy);
    }

    scale_to_stored_basis(ux, uy, mxx, mxy, myy);
    moment_space_collision(ux, uy, mxx, mxy, myy);

    S.global_rho[global_idx] = rho - Geometry::RHO_0;
    S.global_ux[global_idx] = ux;
    S.global_uy[global_idx] = uy;
    S.global_mxx[global_idx] = mxx;
    S.global_mxy[global_idx] = mxy;
    S.global_myy[global_idx] = myy;
}

__global__ void swap_layers(LBMState S, DomainTags T, int y0)
{
    const int x = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
                  static_cast<int>(threadIdx.x);

    if (x >= Geometry::NX)
        return;

    const size_t idx_0 = x + 0 * Geometry::NX;
    const size_t idx_1 = x + 1 * Geometry::NX;
    const size_t idx_2 = x + 2 * Geometry::NX;

    S.d_rho[idx_0] = S.d_rho[idx_1];
    S.d_rho[idx_1] = S.d_rho[idx_2];

    S.d_ux[idx_0] = S.d_ux[idx_1];
    S.d_ux[idx_1] = S.d_ux[idx_2];

    S.d_uy[idx_0] = S.d_uy[idx_1];
    S.d_uy[idx_1] = S.d_uy[idx_2];

    S.d_mxx[idx_0] = S.d_mxx[idx_1];
    S.d_mxx[idx_1] = S.d_mxx[idx_2];

    S.d_mxy[idx_0] = S.d_mxy[idx_1];
    S.d_mxy[idx_1] = S.d_mxy[idx_2];

    S.d_myy[idx_0] = S.d_myy[idx_1];
    S.d_myy[idx_1] = S.d_myy[idx_2];

    if ((y0 + 3) < Geometry::NY)
    {
        const size_t idx_next = x + (y0 + 3) * Geometry::NX;

        S.d_rho[idx_2] = S.global_rho[idx_next];
        S.d_ux[idx_2] = S.global_ux[idx_next];
        S.d_uy[idx_2] = S.global_uy[idx_next];
        S.d_mxx[idx_2] = S.global_mxx[idx_next];
        S.d_mxy[idx_2] = S.global_mxy[idx_next];
        S.d_myy[idx_2] = S.global_myy[idx_next];
    }
}

__global__ void evaluate_last_layer(LBMState S, DomainTags T)
{
    const int x = static_cast<int>(blockIdx.x) * static_cast<int>(blockDim.x) +
                  static_cast<int>(threadIdx.x);

    if (x >= Geometry::NX)
        return;

    int idx = x + (Geometry::NY - 1) * Geometry::NX;

    const uint8_t wall_id = T.d_node[idx];
    const uint32_t valid_ms = T.d_valid[idx];

    real_t pop[Stencil::Q];

    reconstruct_streamed_pop(pop, S, x, 2);

    real_t rho, ux, uy, mxx, mxy, myy;

    if (wall_id != to_u8(NodeId::FLUID))
    {
        Geometry::bc_velocity(x, (Geometry::NY - 1), ux, uy);

        apply_boundary_dirichlet(pop, valid_ms, rho, ux, uy, mxx, mxy, myy);
    }

    scale_to_stored_basis(ux, uy, mxx, mxy, myy);
    moment_space_collision(ux, uy, mxx, mxy, myy);

    S.global_rho[idx] = rho - Geometry::RHO_0;
    S.global_ux[idx] = ux;
    S.global_uy[idx] = uy;
    S.global_mxx[idx] = mxx;
    S.global_mxy[idx] = mxy;
    S.global_myy[idx] = myy;
}

void lbm_mom_step(LBMState &S, const CudaConfig &cfg, const DomainTags &T)
{
    const int bx = (Geometry::NX < 1024) ? Geometry::NX : 1024;
    dim3 block(bx, 1, 1);
    dim3 grid((Geometry::NX + block.x - 1) / block.x, 1, 1);

    load_layers<<<grid, block>>>(S);

    evaluate_first_layer<<<grid, block>>>(S, T);

    for (int y0 = 0; y0 < Geometry::NY - 2; ++y0)
    {
        evaluate_middle_layer<<<grid, block>>>(S, T, y0);

        if (y0 < Geometry::NY - 3)
        {
            swap_layers<<<grid, block>>>(S, T, y0);
        }
    }

    evaluate_last_layer<<<grid, block>>>(S, T);

    CUDA_CHECK(cudaGetLastError());
}
