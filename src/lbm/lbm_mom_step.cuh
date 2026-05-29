#pragma once

#include "state/lbm_state.cuh"
#include "../app/run_context.cuh"
#include "domain/domain_tags.cuh"

#include <vector>

__host__ void lbm_mom_step(std::vector<LBMState> &S,
                           const std::vector<DomainTags> &T,
                           const app::RunContext &ctx);
