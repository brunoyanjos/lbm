#include <iostream>
#include <string>
#include <cuda_runtime.h>

#include "core/cuda_utils.cuh"
#include "app/cuda_config.cuh"
#include "app/simulation.cuh"
#include "app/simulation_summary.cuh"
#include "app/run_context.cuh"

static std::string get_arg(int argc, char **argv, const std::string &key, const std::string &def)
{
    for (int i = 1; i + 1 < argc; ++i)
        if (key == argv[i])
            return argv[i + 1];
    return def;
}

static int get_arg_int(int argc, char **argv, const std::string &key, int def)
{
    const std::string s = get_arg(argc, argv, key, "");
    if (s.empty())
        return def;
    try
    {
        return std::stoi(s);
    }
    catch (...)
    {
        return def;
    }
}

static void configure_simulation_from_checkpoint(app::RunContext &ctx)
{
    if (!ctx.restart_from_checkpoint)
        return;

    std::cout << "[CHECKPOINT] restart requested\n"
              << "[CHECKPOINT] run_id=" << ctx.checkpoint_run_id << "\n"
              << "[CHECKPOINT] dir=" << ctx.checkpoint_dir << "\n"
              << "[CHECKPOINT] metadata/state loading will be wired here.\n";
}

int main(int argc, char **argv)
{
    CUDA_CHECK(cudaSetDevice(0));

    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));

    app::RunContext ctx;
    ctx.out_dir = get_arg(argc, argv, "--out", "runs/default");
    ctx.restart_from_checkpoint = (get_arg_int(argc, argv, "--restart", 0) != 0);
    ctx.checkpoint_run_id = get_arg(argc, argv, "--checkpoint_run_id", "");
    ctx.checkpoint_dir = get_arg(argc, argv, "--checkpoint_dir", "");
    ctx.enable_io = (get_arg_int(argc, argv, "--io", 1) != 0);
    ctx.warmup_steps = get_arg_int(argc, argv, "--warmup", 100);
    ctx.vti_interval = get_arg_int(argc, argv, "--vti_interval", 0);
    ctx.verbose = (get_arg_int(argc, argv, "--verbose", 0) != 0);

    ctx.show_progress = (get_arg_int(argc, argv, "--progress", 1) != 0);
    {
        const int hz = get_arg_int(argc, argv, "--progress_hz", 2);
        ctx.progress_hz = (hz > 0 ? double(hz) : 2.0);
    }

    if (!ctx.enable_io)
        ctx.show_progress = false;

    configure_simulation_from_checkpoint(ctx);

    CudaConfig cfg = make_config();
    if (ctx.restart_from_checkpoint)
    {
        std::cout << "[CHECKPOINT] simulation summary skipped until checkpoint metadata is loaded.\n";
    }
    else
    {
        print_simulation_summary(cfg, prop);
    }

    app::run(cfg, ctx);
    return 0;
}
