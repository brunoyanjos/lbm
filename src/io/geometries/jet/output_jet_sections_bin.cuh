#pragma once
#include <string>

#include "../../../lbm/state/lbm_state.cuh"
#include "../../../lbm/domain/domain_tags.cuh"

namespace io
{
    void write_jet_sections(const LBMState &state, int t, const std::string &out_dir, const DomainTags &tags);
}