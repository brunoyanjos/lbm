#pragma once
#include <cstddef>

constexpr int NX = 1024;
constexpr int NY = NX / 2;

constexpr int JET_Y0 = NY / 2 - 4;
constexpr int JET_Y1 = NY / 2 + 4;

constexpr int JET_H = (JET_Y1 - JET_Y0);