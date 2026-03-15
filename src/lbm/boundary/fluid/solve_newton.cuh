// #pragma once

// #include "../../../core/types.cuh"
// #include "../../../core/math_utils.cuh"
// #include "../../../core/linear_solver.cuh"
// #include "newton_coefficients.cuh"

// __device__ __forceinline__ real_t rel_step(real_t x_new, real_t x_old)
// {
//     const real_t eps = real_t(1e-12);
//     const real_t denom = fmax(r_abs(x_new), eps);
//     return r_abs(x_new - x_old) / denom;
// }

// __device__ __forceinline__ real_t eval_fluid_row(
//     const real_t *__restrict__ row,
//     real_t b_row,
//     real_t ux,
//     real_t uy,
//     real_t mxx,
//     real_t mxy,
//     real_t myy)
// {
//     return row[0] * ux + row[1] * uy +
//            row[2] * ux * ux + row[3] * ux * uy + row[4] * uy * uy +
//            row[5] * mxx + row[6] * mxy + row[7] * myy -
//            b_row;
// }

// __device__ __forceinline__ void build_fluid_jacobian(
//     const FluidNewtonCoefficients &C,
//     real_t ux,
//     real_t uy,
//     real_t *__restrict__ J)
// {
//     auto U = [&](int r, int c) -> real_t &
//     { return J[r * 5 + c]; };

// #pragma unroll
//     for (int i = 0; i < 5; ++i)
//     {
//         U(i, 0) = C.coeff(i, 0) + r::two * C.coeff(i, 2) * ux + C.coeff(i, 3) * uy;
//         U(i, 1) = C.coeff(i, 1) + C.coeff(i, 3) * ux + r::two * C.coeff(i, 4) * uy;
//         U(i, 2) = C.coeff(i, 5);
//         U(i, 3) = C.coeff(i, 6);
//         U(i, 4) = C.coeff(i, 7);
//     }
// }

// __device__ __forceinline__ void build_fluid_rhs(
//     const FluidNewtonCoefficients &C,
//     real_t ux,
//     real_t uy,
//     real_t mxx,
//     real_t mxy,
//     real_t myy,
//     real_t *__restrict__ rhs)
// {
// #pragma unroll
//     for (int i = 0; i < 5; ++i)
//         rhs[i] = -eval_fluid_row(&C.A[i * 8], C.b[i], ux, uy, mxx, mxy, myy);
// }

// __device__ __forceinline__ void solve_fluid_newton(
//     const FluidNewtonCoefficients &C,
//     real_t &ux,
//     real_t &uy,
//     real_t &mxx,
//     real_t &mxy,
//     real_t &myy)
// {
//     real_t J[25]{};
//     real_t rhs[5]{};
//     real_t delta[5]{};

//     real_t error = r::one;
//     int it = 0;
//     constexpr int it_max = 50;
//     constexpr real_t tol = real_t(1e-6);

//     while (error > tol && it++ < it_max)
//     {
//         const real_t ux_old = ux;
//         const real_t uy_old = uy;
//         const real_t mxx_old = mxx;
//         const real_t mxy_old = mxy;
//         const real_t myy_old = myy;

//         build_fluid_jacobian(C, ux, uy, J);
//         build_fluid_rhs(C, ux, uy, mxx, mxy, myy, rhs);

//         // gaussianElimination<5>(J, rhs, delta);

//         ux = ux_old + delta[0];
//         uy = uy_old + delta[1];
//         mxx = mxx_old + delta[2];
//         mxy = mxy_old + delta[3];
//         myy = myy_old + delta[4];

//         error = rel_step(ux, ux_old);
//         error = fmax(error, rel_step(uy, uy_old));
//         error = fmax(error, rel_step(mxx, mxx_old));
//         error = fmax(error, rel_step(mxy, mxy_old));
//         error = fmax(error, rel_step(myy, myy_old));
//     }
// }