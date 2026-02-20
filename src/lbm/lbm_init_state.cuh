#pragma once

#include "state/lbm_state.cuh"
#include "../app/cuda_config.cuh"

__host__ void init_state(LBMState &S, const CudaConfig &cfg);
__host__ void upload_state_to_host(LBMState &S);