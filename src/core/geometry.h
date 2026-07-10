#pragma once
#include <cstddef>

#ifndef LBM_NY
#define LBM_NY 128
#endif

#ifndef LBM_NX
#define LBM_NX (4 * LBM_NY)
#endif

constexpr int NX = LBM_NX;
constexpr int NY = LBM_NY;
