#include <iostream>
#include <string>
#include <cuda_runtime.h>

#include "src/core/cuda_utils.cuh"
#include "src/app/cuda_config.cuh"
#include "src/app/simulation.cuh"
#include "src/app/simulation_summary.cuh"
#include "src/app/run_context.cuh"

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

int main(int argc, char **argv)
{
    CUDA_CHECK(cudaSetDevice(0));

    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));

    CudaConfig cfg = make_config(prop.sharedMemPerBlock);
    print_simulation_summary(cfg, prop);

    app::RunContext ctx;
    ctx.out_dir = get_arg(argc, argv, "--out", "runs/default");
    ctx.enable_io = (get_arg_int(argc, argv, "--io", 1) != 0);
    ctx.warmup_steps = get_arg_int(argc, argv, "--warmup", 100);
    ctx.verbose = (get_arg_int(argc, argv, "--verbose", 0) != 0);

    ctx.show_progress = (get_arg_int(argc, argv, "--progress", 1) != 0);
    {
        const int hz = get_arg_int(argc, argv, "--progress_hz", 2);
        ctx.progress_hz = (hz > 0 ? double(hz) : 2.0);
    }

    if (!ctx.enable_io)
        ctx.show_progress = false;

    app::run(cfg, ctx);
    return 0;
}
