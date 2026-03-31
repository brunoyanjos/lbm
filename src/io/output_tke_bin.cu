#include "output_tke_bin.cuh"

#include "../core/physics.h"
#include "../core/simulation_config.h"
#include "../lbm/stencil_active.cuh"

#include "lbm/hermite/hermite.cuh"
#include "lbm/moment/moment_id.cuh"
#include "core/indexing.cuh"

#include <fstream>
#include <cstdint>
#include <string>

namespace io
{
    static std::string tke_bin_path(const std::string &out_dir)
    {
        return out_dir + "/outputs/tke.bin";
    }

    real_t compute_ke_host_2d(const LBMState &state, const DomainTags &T)
    {
        real_t sum = 0.0;
        int count = 0;

        for (int y = 0; y < NY; ++y)
        {
            for (int x = 0; x < NX; ++x)
            {
                const size_t idx = idxGlobal(x, y);

                if (T.h_node && T.h_node[idx] == to_u8(NodeId::SOLID))
                    continue;

                for (int i = 0; i < Stencil::Q; ++i)
                {
                    const int cx = Stencil::cx(i);
                    const int cy = Stencil::cy(i);

                    const real_t rho = state.h_rho[idx] + RHO_0;
                    const real_t ux = state.h_ux[idx];
                    const real_t uy = state.h_uy[idx];
                    const real_t mxx = state.h_mxx[idx];
                    const real_t mxy = state.h_mxy[idx];
                    const real_t myy = state.h_myy[idx];

                    const real_t fi = Stencil::w(i) * rho *
                                      (r::one +
                                       ux * hermite<MomentId::ux>(i) + uy * hermite<MomentId::uy>(i) +
                                       mxx * hermite<MomentId::mxx>(i) + mxy * hermite<MomentId::mxy>(i) +
                                       myy * hermite<MomentId::myy>(i));

                    const real_t ci2 = r_cast(cx) * r_cast(cx) + r_cast(cy) * r_cast(cy);

                    sum += fi * ci2 * r::half;
                }

                count++;
            }
        }

        real_t norm = count * U_MAX * U_MAX;
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
