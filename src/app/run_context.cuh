#pragma once
#include <string>

namespace app
{

    struct RunContext
    {
        std::string out_dir;   // runs/<id>
        bool enable_io = true; // desliga VTK/upload para benchmark
        int warmup_steps = 100;
        bool verbose = false;

        bool show_progress = true;
        double progress_hz = 2.0;
    };

} // namespace app
