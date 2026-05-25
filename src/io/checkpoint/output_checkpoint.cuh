#pragma once

#include "checkpoint_config.cuh"

#include "../../lbm/state/lbm_state.cuh"

#include <string>

namespace io
{
    void write_checkpoint_current(const LBMState &state, int step, const std::string &out_dir);

    CheckpointConfig read_checkpoint_current(LBMState &state, const std::string &checkpoint_path_or_dir);
}
