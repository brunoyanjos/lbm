#pragma once
#include <string>

#include "../lbm/state/lbm_state.cuh"
#include "../lbm/domain/domain_tags.cuh"

namespace io
{
    void write_annul_profile(const LBMState &state, int t, const std::string &out_dir, const DomainTags &tags);
}