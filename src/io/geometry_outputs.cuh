#pragma once
#include <string>

#include "../lbm/state/lbm_state.cuh"
#include "../lbm/domain/domain_tags.cuh"
#include "../app/cuda_config.cuh"

namespace io
{
    // Chamado dentro do loop (SAVE_INTERVAL)
    void geometry_step_outputs(const LBMState &state,
                               const CudaConfig &cfg,
                               int t,
                               const std::string &out_dir,
                               const DomainTags &tags);

    // Chamado no final (após upload_state_to_host)
    void geometry_final_outputs(const LBMState &state,
                                const CudaConfig &cfg,
                                int t_end,
                                const std::string &out_dir,
                                const DomainTags &tags);
}