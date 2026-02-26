#include "output_tke_bin.cuh"

#include "../core/physics.h"
#include "../core/simulation_config.h"
#include "../lbm/stencil_active.cuh"

#include <fstream>
#include <cstdint>
#include <string>

namespace io
{
    static std::string tke_bin_path(const std::string &out_dir)
    {
        return out_dir + "/outputs/tke.bin";
    }

    real_t compute_ke_host_2d(const LBMState &state)
    {
        real_t sum = 0.0;

        for (int i = 0; i < int(state.N); ++i)
        {
            const real_t ux = (real_t)state.h_ux[i] / Stencil::as2;
            const real_t uy = (real_t)state.h_uy[i] / Stencil::as2;
            sum += 0.5 * (ux * ux + uy * uy);
        }

        real_t norm = state.N * U_LID * U_LID;
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
}
