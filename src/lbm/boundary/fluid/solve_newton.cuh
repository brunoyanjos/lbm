#pragma once

#include "core/math_utils.cuh"
#include "lbm/moment/node_moments.cuh"

#include "lbm/boundary/common/system_data.cuh"

namespace boundary::fluid
{
    __device__ __forceinline__ real_t rel_step(real_t x_new, real_t x_old)
    {
        const real_t eps = real_t(1e-12);
        const real_t denom = fmax(r_abs(x_new), eps);
        return r_abs(x_new - x_old) / denom;
    }

    __device__ __forceinline__
        real_t
        max_rel_step(const NodeMoments &a, const NodeMoments &b)
    {
        real_t error = rel_step(a.ux, b.ux);
        error = fmax(error, rel_step(a.uy, b.uy));
        error = fmax(error, rel_step(a.mxx, b.mxx));
        error = fmax(error, rel_step(a.mxy, b.mxy));
        error = fmax(error, rel_step(a.myy, b.myy));
        return error;
    }

    template <std::size_t N, std::size_t E, int RegOrder, bool Rec, bool HighOrder>
    __device__ __forceinline__ real_t eval_row(const SystemData<N, E> &S,
                                               const NodeMomentsFor<RegOrder, Rec, HighOrder> &M,
                                               int i)
    {
        const real_t ux2 = M.ux * M.ux;
        const real_t uxuy = M.ux * M.uy;
        const real_t uy2 = M.uy * M.uy;

        real_t value = M.ux * S.coeff(i, 0) + M.uy * S.coeff(i, 1) +
                       ux2 * S.coeff(i, 2) + uxuy * S.coeff(i, 3) + uy2 * S.coeff(i, 4) +
                       M.mxx * S.coeff(i, 5) + M.mxy * S.coeff(i, 6) + M.myy * S.coeff(i, 7) -
                       S.b[i];

        if constexpr (RegOrder >= 3)
        {
            const real_t ux3 = M.ux * M.ux * M.ux;
            const real_t ux2uy = M.ux * M.ux * M.uy;
            const real_t uxuy2 = M.ux * M.uy * M.uy;
            const real_t uy3 = M.uy * M.uy * M.uy;

            value += ux3 * S.coeff(i, 8) + ux2uy * S.coeff(i, 9) +
                     uxuy2 * S.coeff(i, 10) + uy3 * S.coeff(i, 11);

            if constexpr (Rec)
            {
                value += S.coeff(i, 12) * M.ux * M.mxx + S.coeff(i, 13) * M.uy * M.mxx +
                         S.coeff(i, 14) * M.ux * M.mxy + S.coeff(i, 15) * M.uy * M.mxy +
                         S.coeff(i, 16) * M.ux * M.myy + S.coeff(i, 17) * M.uy * M.myy;
            }
        }

        return value;
    }

    template <std::size_t N, std::size_t E, int RegOrder, bool Rec, bool HighOrder>
    __device__ __forceinline__ void build_newton_step(SystemData<N> &G,
                                                      const SystemData<N, E> &S,
                                                      const NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
    {
#pragma unroll
        for (int i = 0; i < int(N); ++i)
        {
            G.coeff(i, 0) = S.coeff(i, 0) + r::two * S.coeff(i, 2) * M.ux + S.coeff(i, 3) * M.uy;
            G.coeff(i, 1) = S.coeff(i, 1) + S.coeff(i, 3) * M.ux + r::two * S.coeff(i, 4) * M.uy;
            G.coeff(i, 2) = S.coeff(i, 5);
            G.coeff(i, 3) = S.coeff(i, 6);
            G.coeff(i, 4) = S.coeff(i, 7);

            if constexpr (RegOrder >= 3)
            {
                G.coeff(i, 0) += r::three * S.coeff(i, 8) * M.ux * M.ux +
                                 r::two * S.coeff(i, 9) * M.ux * M.uy +
                                 S.coeff(i, 10) * M.uy * M.uy;

                G.coeff(i, 1) += S.coeff(i, 9) * M.ux * M.ux +
                                 r::two * S.coeff(i, 10) * M.ux * M.uy +
                                 r::three * S.coeff(i, 11) * M.uy * M.uy;

                if constexpr (Rec)
                {
                    G.coeff(i, 0) += S.coeff(i, 12) * M.mxx + S.coeff(i, 14) * M.mxy + S.coeff(i, 16) * M.myy;
                    G.coeff(i, 1) += S.coeff(i, 13) * M.mxx + S.coeff(i, 15) * M.mxy + S.coeff(i, 17) * M.myy;
                    G.coeff(i, 2) += S.coeff(i, 12) * M.ux + S.coeff(i, 13) * M.uy;
                    G.coeff(i, 3) += S.coeff(i, 14) * M.ux + S.coeff(i, 15) * M.uy;
                    G.coeff(i, 4) += S.coeff(i, 16) * M.ux + S.coeff(i, 17) * M.uy;
                }
            }

            G.b[i] = -eval_row(S, M, i);
        }
    }

    template <std::size_t N, std::size_t E, int RegOrder, bool Rec, bool HighOrder>
    __device__ __forceinline__ void solve_newton(const SystemData<N, E> &S,
                                                 NodeMomentsFor<RegOrder, Rec, HighOrder> &M)
    {
        SystemData<N> G{};

        real_t error = r::one;
        int it = 0;
        constexpr int it_max = 50;
        constexpr real_t tol = real_t(1e-6);

        while (error > tol && it++ < it_max)
        {
            const NodeMoments old = M;

            build_newton_step(G, S, M);
            gaussianElimination(G);

            M.ux += G.x[0];
            M.uy += G.x[1];
            M.mxx += G.x[2];
            M.mxy += G.x[3];
            M.myy += G.x[4];

            error = max_rel_step(M, old);
        }
    }
}