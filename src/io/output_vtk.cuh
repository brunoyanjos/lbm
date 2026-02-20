#pragma once
#include "../lbm/state/lbm_state.cuh"
#include "../core/geometry.h"
#include <string>

namespace io
{
    __host__ void write_vtk(const LBMState &S, const CudaConfig &cfg, int step, const std::string &out_dir);
}
