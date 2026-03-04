#include "output_jet_sections_bin.cuh"

#if defined(LBM_GEOM_JET)

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include "../../../core/types.cuh"
#include "../../../lbm/stencil_active.cuh"
#include "../../../core/geometries/jet/geometry.h"
#include "../../../core/geometries/jet/physics.h"

namespace io
{
    static std::string jet_sections_path(const std::string &out_dir)
    {
        return out_dir + "/outputs/jet_sections.bin";
    }

    void write_jet_sections(const LBMState &state, int t, const std::string &out_dir, const DomainTags &tags)
    {
        (void)tags;

        const int x_sec[3] = {NX / 4, NX / 2, (3 * NX) / 4};
        const int nsec = 3;

        std::vector<float> ux_all;
        ux_all.resize((size_t)nsec * (size_t)NY);

        for (int k = 0; k < nsec; ++k)
        {
            const int xs = x_sec[k];

            for (int y = 0; y < NY; ++y)
            {
                const int idx = xs + NX * y;
                const int idx_n = (xs - 1) + NX * y;
                const int idx_n_safe = (xs == 0) ? idx : idx_n;

                const real_t ux = (state.h_ux[idx] / Stencil::as2 + state.h_ux[idx_n_safe] / Stencil::as2) * real_t(0.5);
                ux_all[(size_t)k * (size_t)NY + (size_t)y] = (float)ux;
            }
        }

        const std::string path = jet_sections_path(out_dir);
        std::ofstream f(path, std::ios::binary);

        const char magic[4] = {'J', 'S', 'C', '1'};
        f.write(magic, 4);

        const int32_t nx32 = NX;
        const int32_t ny32 = NY;
        const int32_t t32 = (int32_t)t;
        const int32_t ns32 = (int32_t)nsec;

        f.write(reinterpret_cast<const char *>(&nx32), sizeof(nx32));
        f.write(reinterpret_cast<const char *>(&ny32), sizeof(ny32));
        f.write(reinterpret_cast<const char *>(&t32), sizeof(t32));
        f.write(reinterpret_cast<const char *>(&ns32), sizeof(ns32));

        for (int k = 0; k < nsec; ++k)
        {
            const int32_t xs32 = (int32_t)x_sec[k];
            f.write(reinterpret_cast<const char *>(&xs32), sizeof(xs32));
        }

        f.write(reinterpret_cast<const char *>(ux_all.data()), sizeof(float) * ux_all.size());
    }
}

#endif // LBM_GEOM_JET