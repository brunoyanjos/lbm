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
        real_t sum = 0.0;
        int count = 0;

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
}
