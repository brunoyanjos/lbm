#pragma once
#include <string>

#include "../lbm/state/lbm_state.cuh"

namespace io
{
    void write_centerline_profiles(const LBMState &state, int t, const std::string &out_dir);
}