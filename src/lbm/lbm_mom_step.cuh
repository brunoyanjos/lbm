#pragma once

#include "state/lbm_state.cuh"
#include "../app/cuda_config.cuh"
#include "domain/domain_tags.cuh"

__host__ void lbm_mom_step(LBMState &S, const CudaConfig &cfg, const DomainTags &T);