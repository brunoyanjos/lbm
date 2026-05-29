#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cuda_runtime.h>

#include "core/cuda_utils.cuh"
#include "app/cuda_config.cuh"
#include "app/simulation.cuh"
#include "app/simulation_summary.cuh"
#include "app/run_context.cuh"
#include "app/domain_partition.cuh"
#include "io/checkpoint/checkpoint_config.cuh"

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

static std::vector<int> parse_device_list(const std::string &csv)
{
    std::vector<int> devices;
    std::stringstream ss(csv);
    std::string item;

    while (std::getline(ss, item, ','))
    {
        if (item.empty())
            throw std::invalid_argument("empty device id");

        size_t consumed = 0;
        const int device_id = std::stoi(item, &consumed);
        if (consumed != item.size())
            throw std::invalid_argument("invalid device id: " + item);

        devices.push_back(device_id);
    }

    if (devices.empty())
        throw std::invalid_argument("device list is empty");

    return devices;
}

static void configure_simulation_from_checkpoint(app::RunContext &ctx)
{
    if (!ctx.restart_from_checkpoint)
        return;

    std::cout << "[CHECKPOINT] restart requested\n"
              << "[CHECKPOINT] run_id=" << ctx.checkpoint_run_id << "\n"
              << "[CHECKPOINT] dir=" << ctx.checkpoint_dir << "\n"
              << "[CHECKPOINT] state will be loaded before time stepping.\n";
}

int main(int argc, char **argv)
{
    const int device_id = get_arg_int(argc, argv, "--device", 0);
    std::vector<int> device_ids;
    try
    {
        const std::string devices_arg = get_arg(argc, argv, "--devices", "");
        device_ids = devices_arg.empty() ? std::vector<int>{device_id} : parse_device_list(devices_arg);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Invalid --devices value: " << e.what() << "\n";
        return 1;
    }
    const int primary_device_id = device_ids.front();

    int device_count = 0;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    if (device_count <= 0)
    {
        std::cerr << "No CUDA devices found.\n";
        return 1;
    }
    if (static_cast<int>(device_ids.size()) > NY)
    {
        std::cerr << "Too many partitions for NY=" << NY
                  << ": requested " << device_ids.size() << " devices.\n";
        return 1;
    }

    for (int id : device_ids)
    {
        if (id < 0 || id >= device_count)
        {
            std::cerr << "Invalid CUDA device " << id
                      << ". Available device ids: 0.." << (device_count - 1) << "\n";
            return 1;
        }
    }

    CUDA_CHECK(cudaSetDevice(primary_device_id));

    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, primary_device_id));

    app::RunContext ctx;
    ctx.out_dir = get_arg(argc, argv, "--out", "runs/default");
    ctx.restart_from_checkpoint = (get_arg_int(argc, argv, "--restart", 0) != 0);
    ctx.checkpoint_run_id = get_arg(argc, argv, "--checkpoint_run_id", "");
    ctx.checkpoint_dir = get_arg(argc, argv, "--checkpoint_dir", "");
    ctx.enable_io = (get_arg_int(argc, argv, "--io", 1) != 0);
    ctx.warmup_steps = get_arg_int(argc, argv, "--warmup", 100);
    ctx.vti_interval = get_arg_int(argc, argv, "--vti_interval", 0);
    ctx.verbose = (get_arg_int(argc, argv, "--verbose", 0) != 0);
    ctx.partitions = app::make_domain_partitions(device_ids);

    ctx.show_progress = (get_arg_int(argc, argv, "--progress", 1) != 0);
    {
        const int hz = get_arg_int(argc, argv, "--progress_hz", 2);
        ctx.progress_hz = (hz > 0 ? double(hz) : 2.0);
    }

    if (!ctx.enable_io)
        ctx.show_progress = false;

    configure_simulation_from_checkpoint(ctx);

    std::cout << "CUDA primary device: " << primary_device_id << " / " << device_count << "\n";
    if (ctx.partitions.size() > 1)
        std::cout << "Multi-device partitions are planned; kernels still execute on the primary device.\n";
    print_simulation_summary(prop, ctx.partitions);

    app::run(ctx);
    return 0;
}
