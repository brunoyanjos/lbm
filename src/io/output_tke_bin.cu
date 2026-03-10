#include "output_tke_bin.cuh"

#include "../geometries/active_geometry.cuh"
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

    real_t compute_ke_host_2d(const LBMState &state, const uint8_t *h_node)
    {
        real_t sum = real_t(0);
        size_t n_active = 0;

        const uint8_t SOLID = to_u8(NodeId::SOLID);

        for (size_t i = 0; i < state.N; ++i)
        {
            if (h_node && h_node[i] == SOLID)
                continue;

            const real_t ux = (real_t)state.h_ux[i] / Stencil::as2;
            const real_t uy = (real_t)state.h_uy[i] / Stencil::as2;

            sum += real_t(0.5) * (ux * ux + uy * uy);
            ++n_active;
        }

        if (n_active == 0)
            return real_t(0);

        // normalização: média de KE adimensionalizada por U_MAX^2
        const real_t inv_norm = real_t(1) / (real_t(n_active) * Geometry::U_MAX * Geometry::U_MAX);
        return sum * inv_norm;
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
