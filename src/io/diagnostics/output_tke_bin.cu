#include "output_tke_bin.cuh"

#include "../../core/physics.h"
#include "../../core/simulation_config.h"
#include "../../lbm/stencil_active.cuh"
#include "lbm/hermite/hermite.cuh"
#include "lbm/moment/scale_factor.cuh"
#include "core/indexing.cuh"

#include <fstream>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>

namespace io
{
    namespace
    {
        real_t compute_ke_host_2d_impl(const LBMState &state,
                                       FlowAverages *averages)
        {
            real_t sum = 0.0;
            int count = 0;

            const bool sample_averages = (averages != nullptr);
            const double next_sample = sample_averages
                                           ? static_cast<double>(averages->samples + 1)
                                           : 0.0;
            if (sample_averages)
                validate_flow_averages_shape(*averages);

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

                    const double uxd = static_cast<double>(ux) * inv_scale_factor<MomentId::ux>();
                    const double uyd = static_cast<double>(uy) * inv_scale_factor<MomentId::uy>();

                    if (sample_averages)
                    {
                        sample_flow_averages_node(*averages, idx, uxd, uyd, next_sample);
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
                        const real_t cix_ux = r_cast(cx) - static_cast<real_t>(uxd);
                        const real_t ciy_uy = r_cast(cy) - static_cast<real_t>(uyd);

                        const real_t ci_u2 = cix_ux * cix_ux + ciy_uy * ciy_uy;

                        rho_E += fi * ci2 * r::half;
                        rho_e += fi * ci_u2 * r::half;
                    }

                    sum += (rho_E - rho_e) * inv_rho;
                    count++;
                }
            }

            real_t norm = count * U * U;
            real_t inv_norm = real_t(1) / norm;

            sum *= inv_norm;

            if (sample_averages)
                ++averages->samples;

            return sum;
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

    real_t compute_ke_host_2d(const LBMState &state)
    {
        return compute_ke_host_2d_impl(state, nullptr);
    }

    real_t compute_ke_and_sample_flow_averages_host_2d(const LBMState &state,
                                                       FlowAverages &averages)
    {
        return compute_ke_host_2d_impl(state, &averages);
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
}
