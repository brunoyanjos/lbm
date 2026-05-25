#pragma once

#include "../../lbm/state/lbm_state.cuh"

#include <string>

namespace io
{
    void write_checkpoint_current(const LBMState &state, int step, const std::string &out_dir);
}
