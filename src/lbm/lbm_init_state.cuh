#pragma once

#include "state/lbm_state.cuh"
#include "app/run_context.cuh"

#include <vector>

__host__ void init_state(std::vector<LBMState> &S, const app::RunContext &ctx);
__host__ void upload_state_to_host(std::vector<LBMState> &S, const app::RunContext &ctx);
__host__ void upload_state_to_host(LBMState &S);
