#pragma once

#include "../lbm/state/lbm_state.cuh"
#include "../lbm/domain/domain_tags.cuh"
#include "../core/types.cuh"

#include <string>
#include <cstdint>

namespace io
{
    real_t compute_ke_host_2d(const LBMState &state, const uint8_t *h_node);

    void tke_bin_append(const std::string &out_dir, int t, double ke);
}