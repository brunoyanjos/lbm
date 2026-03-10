#pragma once

#include <string>

#include "../../lbm/state/lbm_state.cuh"
#include "../../lbm/domain/domain_tags.cuh"

namespace ANNUL
{
    void outputs(const LBMState &state, int t, const std::string &out_dir, const DomainTags &tags);
}