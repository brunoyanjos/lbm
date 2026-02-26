#pragma once

#include "../lbm/state/lbm_state.cuh"
#include "../core/types.cuh"

#include <string>
#include <cstdint>

namespace io
{
    real_t compute_ke_host_2d(const LBMState &state);

    void tke_bin_append(const std::string &out_dir, int t, double ke);
}