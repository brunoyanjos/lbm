#pragma once

#include <cstdint>
#include <type_traits>

#include "../lbm/stencil_active.cuh"

static_assert(Stencil::Q <= 64, "Stencil::Q must be <= 64 for mask_t");

using mask_t = std::conditional_t<(Stencil::Q <= 32), uint32_t, uint64_t>;