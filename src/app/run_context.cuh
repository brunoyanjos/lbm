#pragma once
#include <string>

namespace app
{

    struct RunContext
    {
        std::string out_dir;
        bool restart_from_checkpoint = false;
        std::string checkpoint_run_id;
        std::string checkpoint_dir;
        bool enable_io = true;
        int warmup_steps = 100;
        int vti_interval = 0;
        bool verbose = false;

        bool show_progress = true;
        double progress_hz = 2.0;
    };

} // namespace app
