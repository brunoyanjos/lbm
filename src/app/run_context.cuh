#pragma once
#include <string>

namespace app
{

    struct RunContext
    {
        std::string out_dir;
        bool enable_io = true;
        int warmup_steps = 100;
        bool verbose = false;

        bool show_progress = true;
        double progress_hz = 2.0;
    };

} // namespace app
