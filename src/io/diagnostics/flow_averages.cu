#include "flow_averages.cuh"

#include "../../core/geometry.h"
#include "../../core/simulation_config.h"

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace io
{
    namespace
    {
        constexpr char AVG_MAGIC[8] = {'L', 'B', 'M', 'A', 'V', 'G', '1', '\0'};

        struct FlowAveragesHeader
        {
            char magic[8];
            std::int32_t nx;
            std::int32_t ny;
            std::int32_t field_count;
            std::int32_t value_bytes;
            std::int64_t node_count;
            std::int64_t samples;
            std::int64_t step;
            std::int64_t t_star;
        };

        static_assert(std::is_trivially_copyable_v<FlowAveragesHeader>,
                      "FlowAveragesHeader must stay binary-serializable");

        std::filesystem::path checkpoint_average_path(const std::string &out_dir, int step)
        {
            namespace fs = std::filesystem;

            const int t_star = step / SAVE_INTERVAL;
            std::ostringstream filename;
            filename << "averages_tstar_" << std::setw(7) << std::setfill('0') << t_star
                     << "_step_" << std::setw(9) << std::setfill('0') << step << ".bin";

            return fs::path(out_dir) / "checkpoints" / filename.str();
        }

        std::filesystem::path source_average_path_from_checkpoint(const std::string &checkpoint_path_or_dir,
                                                                  int checkpoint_step)
        {
            namespace fs = std::filesystem;

            fs::path path(checkpoint_path_or_dir);
            fs::path dir = fs::is_regular_file(path) ? path.parent_path() : path;

            if (dir.filename() != "checkpoints")
            {
                if (fs::exists(dir / "checkpoints"))
                    dir = dir / "checkpoints";
                else
                    dir = dir.parent_path() / "checkpoints";
            }

            return checkpoint_average_path(dir.parent_path().string(), checkpoint_step);
        }

        void write_flow_averages_file(const FlowAverages &averages, int step, const std::filesystem::path &path)
        {
            namespace fs = std::filesystem;
            validate_flow_averages_shape(averages);
            fs::create_directories(path.parent_path());

            std::ofstream f(path, std::ios::binary);
            if (!f.is_open())
            {
                std::cerr << "Could not open flow averages file for writing: " << path.string() << "\n";
                return;
            }

            FlowAveragesHeader header{};
            std::memcpy(header.magic, AVG_MAGIC, sizeof(header.magic));
            header.nx = NX;
            header.ny = NY;
            header.field_count = 5;
            header.value_bytes = sizeof(double);
            header.node_count = static_cast<std::int64_t>(NX) * static_cast<std::int64_t>(NY);
            header.samples = averages.samples;
            header.step = step;
            header.t_star = step / SAVE_INTERVAL;

            f.write(reinterpret_cast<const char *>(&header), sizeof(header));
            f.write(reinterpret_cast<const char *>(averages.ux.data()),
                    static_cast<std::streamsize>(averages.ux.size() * sizeof(double)));
            f.write(reinterpret_cast<const char *>(averages.uy.data()),
                    static_cast<std::streamsize>(averages.uy.size() * sizeof(double)));
            f.write(reinterpret_cast<const char *>(averages.uxux.data()),
                    static_cast<std::streamsize>(averages.uxux.size() * sizeof(double)));
            f.write(reinterpret_cast<const char *>(averages.uxuy.data()),
                    static_cast<std::streamsize>(averages.uxuy.size() * sizeof(double)));
            f.write(reinterpret_cast<const char *>(averages.uyuy.data()),
                    static_cast<std::streamsize>(averages.uyuy.size() * sizeof(double)));

            if (!f.good())
                std::cerr << "Flow averages write failed: " << path.string() << "\n";
        }
    }

    FlowAverages make_flow_averages(std::size_t node_count)
    {
        FlowAverages averages{};
        averages.ux.assign(node_count, 0.0);
        averages.uy.assign(node_count, 0.0);
        averages.uxux.assign(node_count, 0.0);
        averages.uxuy.assign(node_count, 0.0);
        averages.uyuy.assign(node_count, 0.0);
        return averages;
    }

    void validate_flow_averages_shape(const FlowAverages &averages)
    {
        const std::size_t node_count = static_cast<std::size_t>(NX) * static_cast<std::size_t>(NY);
        if (averages.ux.size() != node_count ||
            averages.uy.size() != node_count ||
            averages.uxux.size() != node_count ||
            averages.uxuy.size() != node_count ||
            averages.uyuy.size() != node_count)
        {
            throw std::runtime_error("FlowAverages field sizes do not match the active grid");
        }
    }

    void sample_flow_averages_node(FlowAverages &averages,
                                   std::size_t idx,
                                   double ux,
                                   double uy,
                                   double next_sample)
    {
        averages.ux[idx] += (ux - averages.ux[idx]) / next_sample;
        averages.uy[idx] += (uy - averages.uy[idx]) / next_sample;
        averages.uxux[idx] += (ux * ux - averages.uxux[idx]) / next_sample;
        averages.uxuy[idx] += (ux * uy - averages.uxuy[idx]) / next_sample;
        averages.uyuy[idx] += (uy * uy - averages.uyuy[idx]) / next_sample;
    }

    bool seed_flow_averages_from_checkpoint(const std::string &checkpoint_path_or_dir,
                                            FlowAverages &averages,
                                            int checkpoint_step)
    {
        namespace fs = std::filesystem;

        const fs::path src_path = source_average_path_from_checkpoint(checkpoint_path_or_dir,
                                                                      checkpoint_step);
        if (!fs::exists(src_path))
        {
            std::cerr << "[CHECKPOINT] previous flow averages not found: " << src_path.string() << "\n";
            return false;
        }

        std::ifstream f(src_path, std::ios::binary);
        if (!f.is_open())
        {
            std::cerr << "[CHECKPOINT] could not seed flow averages from: " << src_path.string() << "\n";
            return false;
        }

        FlowAveragesHeader header{};
        f.read(reinterpret_cast<char *>(&header), sizeof(header));
        if (!f.good() || std::memcmp(header.magic, AVG_MAGIC, sizeof(header.magic)) != 0)
        {
            std::cerr << "[CHECKPOINT] invalid flow averages file: " << src_path.string() << "\n";
            return false;
        }

        const std::int64_t node_count = static_cast<std::int64_t>(NX) * static_cast<std::int64_t>(NY);
        if (header.nx != NX || header.ny != NY ||
            header.field_count != 5 ||
            header.value_bytes != static_cast<std::int32_t>(sizeof(double)) ||
            header.node_count != node_count)
        {
            std::cerr << "[CHECKPOINT] flow averages shape does not match this executable: "
                      << src_path.string() << "\n";
            return false;
        }

        averages = make_flow_averages(static_cast<std::size_t>(node_count));
        averages.samples = header.samples;

        f.read(reinterpret_cast<char *>(averages.ux.data()),
               static_cast<std::streamsize>(averages.ux.size() * sizeof(double)));
        f.read(reinterpret_cast<char *>(averages.uy.data()),
               static_cast<std::streamsize>(averages.uy.size() * sizeof(double)));
        f.read(reinterpret_cast<char *>(averages.uxux.data()),
               static_cast<std::streamsize>(averages.uxux.size() * sizeof(double)));
        f.read(reinterpret_cast<char *>(averages.uxuy.data()),
               static_cast<std::streamsize>(averages.uxuy.size() * sizeof(double)));
        f.read(reinterpret_cast<char *>(averages.uyuy.data()),
               static_cast<std::streamsize>(averages.uyuy.size() * sizeof(double)));

        if (!f.good())
        {
            std::cerr << "[CHECKPOINT] could not read flow averages payload: " << src_path.string() << "\n";
            averages = make_flow_averages(static_cast<std::size_t>(node_count));
            return false;
        }

        std::cout << "[CHECKPOINT] seeded flow averages with " << averages.samples
                  << " samples from " << src_path.string() << "\n";
        return true;
    }

    void write_flow_averages(const FlowAverages &averages,
                             int step,
                             const std::string &out_dir)
    {
        write_flow_averages_file(averages, step,
                                 std::filesystem::path(out_dir) / "outputs" / "flow_averages.bin");
    }

    void write_flow_averages_checkpoint(const FlowAverages &averages,
                                        int step,
                                        const std::string &out_dir)
    {
        write_flow_averages_file(averages, step, checkpoint_average_path(out_dir, step));
    }
}
