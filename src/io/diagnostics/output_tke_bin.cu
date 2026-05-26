#include "output_tke_bin.cuh"

#include "../../core/physics.h"
#include "../../core/simulation_config.h"
#include "../../lbm/stencil_active.cuh"
#include "lbm/hermite/hermite.cuh"
#include "lbm/moment/scale_factor.cuh"
#include "core/indexing.cuh"

#include <fstream>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
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

    static std::string tke_bin_path(const std::string &out_dir)
    {
        return out_dir + "/outputs/tke.bin";
    }

    static std::filesystem::path source_tke_path_from_checkpoint(const std::string &checkpoint_path_or_dir)
    {
        namespace fs = std::filesystem;

        fs::path path(checkpoint_path_or_dir);
        fs::path dir = fs::is_regular_file(path) ? path.parent_path() : path;

        if (dir.filename() == "checkpoints")
            return dir.parent_path() / "outputs" / "tke.bin";

        if (fs::exists(dir / "outputs" / "tke.bin"))
            return dir / "outputs" / "tke.bin";

        return dir.parent_path() / "outputs" / "tke.bin";
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

    real_t compute_ke_host_2d(const LBMState &state)
    {
        FlowAverages unused = make_flow_averages(0);
        return compute_ke_and_accumulate_flow_averages_host_2d(state, unused);
    }

    real_t compute_ke_and_accumulate_flow_averages_host_2d(const LBMState &state,
                                                           FlowAverages &averages)
    {
        real_t sum = 0.0;
        int count = 0;
        const bool accumulate_averages = !averages.ux.empty();
        const double next_sample = static_cast<double>(averages.samples + 1);
        if (accumulate_averages)
            validate_flow_averages_shape(averages);

        for (int y = 0; y < NY; ++y)
        {
            for (int x = 0; x < NX; ++x)
            {
                const size_t idx = idxGlobal(x, y);

                real_t rho_E = 0.0;
                real_t rho_e = 0.0;

                const real_t rho = state.h_rho[idx] + RHO_0;
                const real_t inv_rho = r::one / rho;
                const real_t ux = state.h_ux[idx];
                const real_t uy = state.h_uy[idx];
                const real_t mxx = state.h_mxx[idx];
                const real_t mxy = state.h_mxy[idx];
                const real_t myy = state.h_myy[idx];

                if (accumulate_averages)
                {
                    const double uxd = static_cast<double>(ux) * inv_scale_factor<MomentId::ux>();
                    const double uyd = static_cast<double>(uy) * inv_scale_factor<MomentId::uy>();
                    averages.ux[idx] += (uxd - averages.ux[idx]) / next_sample;
                    averages.uy[idx] += (uyd - averages.uy[idx]) / next_sample;
                    averages.uxux[idx] += (uxd * uxd - averages.uxux[idx]) / next_sample;
                    averages.uxuy[idx] += (uxd * uyd - averages.uxuy[idx]) / next_sample;
                    averages.uyuy[idx] += (uyd * uyd - averages.uyuy[idx]) / next_sample;
                }

                for (int i = 0; i < Stencil::Q; ++i)
                {
                    const int cx = Stencil::cx(i);
                    const int cy = Stencil::cy(i);

                    const real_t fi = Stencil::w(i) * rho *
                                      (r::one +
                                       ux * hermite<MomentId::ux>(i) + uy * hermite<MomentId::uy>(i) +
                                       mxx * hermite<MomentId::mxx>(i) + mxy * hermite<MomentId::mxy>(i) +
                                       myy * hermite<MomentId::myy>(i));

                    const real_t ci2 = r_cast(cx) * r_cast(cx) + r_cast(cy) * r_cast(cy);
                    const real_t cix_ux = r_cast(cx) - ux * inv_scale_factor<MomentId::ux>();
                    const real_t ciy_uy = r_cast(cy) - uy * inv_scale_factor<MomentId::uy>();

                    const real_t ci_u2 = cix_ux * cix_ux + ciy_uy * ciy_uy;

                    rho_E += fi * ci2 * r::half;
                    rho_e += fi * ci_u2 * r::half;
                }

                sum += (rho_E - rho_e) * inv_rho;

                count++;
            }
        }

        real_t norm = count * U_LID * U_LID;
        real_t inv_norm = real_t(1) / norm;

        sum *= inv_norm;

        if (accumulate_averages)
            ++averages.samples;

        return sum;
    }

    void tke_bin_append(const std::string &out_dir, int t, double ke)
    {
        const std::string path = tke_bin_path(out_dir);
        std::ofstream f(path, std::ios::binary | std::ios::app);

        const int64_t tstar = int64_t(t) / int64_t(SAVE_INTERVAL);

        f.write(reinterpret_cast<const char *>(&tstar), sizeof(tstar));
        f.write(reinterpret_cast<const char *>(&ke), sizeof(ke));
    }

    void seed_tke_history_from_checkpoint(const std::string &checkpoint_path_or_dir,
                                          const std::string &out_dir,
                                          int checkpoint_step)
    {
        namespace fs = std::filesystem;

        const fs::path src_path = source_tke_path_from_checkpoint(checkpoint_path_or_dir);
        if (!fs::exists(src_path))
        {
            std::cerr << "[CHECKPOINT] previous TKE history not found: " << src_path.string() << "\n";
            return;
        }

        const fs::path dst_path = fs::path(out_dir) / "outputs" / "tke.bin";
        fs::create_directories(dst_path.parent_path());

        if (fs::exists(dst_path) && fs::equivalent(src_path, dst_path))
            return;

        std::ifstream src(src_path, std::ios::binary);
        std::ofstream dst(dst_path, std::ios::binary | std::ios::trunc);
        if (!src.is_open() || !dst.is_open())
        {
            std::cerr << "[CHECKPOINT] could not seed TKE history from: " << src_path.string() << "\n";
            return;
        }

        const int64_t checkpoint_tstar = int64_t(checkpoint_step) / int64_t(SAVE_INTERVAL);
        int64_t tstar = 0;
        double ke = 0.0;
        int64_t copied = 0;

        while (src.read(reinterpret_cast<char *>(&tstar), sizeof(tstar)))
        {
            if (!src.read(reinterpret_cast<char *>(&ke), sizeof(ke)))
                break;

            if (tstar > checkpoint_tstar)
                break;

            dst.write(reinterpret_cast<const char *>(&tstar), sizeof(tstar));
            dst.write(reinterpret_cast<const char *>(&ke), sizeof(ke));
            ++copied;
        }

        std::cout << "[CHECKPOINT] seeded " << copied
                  << " TKE records from " << src_path.string() << "\n";
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
