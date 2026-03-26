#pragma once
#include "../lbm/state/lbm_state.cuh"
#include "../core/geometry.h"
#include "../lbm/domain/domain_tags.cuh"
#include <string>

namespace io
{
    __host__ void write_vti(const LBMState &S, int step, const std::string &out_dir);
}
