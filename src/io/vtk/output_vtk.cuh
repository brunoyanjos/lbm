#pragma once
#include "../../lbm/state/lbm_state.cuh"
#include "../../core/geometry.h"
#include <string>
#include <vector>

namespace io
{
    __host__ void write_vti(const std::vector<LBMState> &states, int step, const std::string &out_dir);
}
