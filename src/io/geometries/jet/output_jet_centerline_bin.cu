#include "output_jet_centerline_bin.cuh"

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
    static std::string jet_centerline_path(const std::string &out_dir)
    {
        return out_dir + "/outputs/jet_centerline.bin";
    }

    void write_jet_centerline(const LBMState &state, int t, const std::string &out_dir, const DomainTags &tags)
    {
        (void)tags;

        // se existir JET_Y0/JET_Y1, use o centro da janela; senão, meio do domínio
#ifdef JET_Y0
        const int y_line = (JET_Y0 + JET_Y1) / 2;
#else
        const int y_line = NY / 2;
#endif

        std::vector<float> ux_x(NX);

        for (int x = 0; x < NX; ++x)
        {
            const int idx = x + NX * y_line;
            const int idx_n = (x - 1) + NX * y_line;
            const int idx_n_safe = (x == 0) ? idx : idx_n;

            const real_t ux = (state.h_ux[idx] / Stencil::as2 + state.h_ux[idx_n_safe] / Stencil::as2) * real_t(0.5);
            ux_x[x] = (float)ux;
        }

        const std::string path = jet_centerline_path(out_dir);
        std::ofstream f(path, std::ios::binary);

        const char magic[4] = {'J', 'C', 'L', '1'};
        f.write(magic, 4);

        const int32_t nx32 = NX;
        const int32_t ny32 = NY;
        const int32_t t32 = (int32_t)t;
        const int32_t yl32 = (int32_t)y_line;

        f.write(reinterpret_cast<const char *>(&nx32), sizeof(nx32));
        f.write(reinterpret_cast<const char *>(&ny32), sizeof(ny32));
        f.write(reinterpret_cast<const char *>(&t32), sizeof(t32));
        f.write(reinterpret_cast<const char *>(&yl32), sizeof(yl32));

        f.write(reinterpret_cast<const char *>(ux_x.data()), sizeof(float) * ux_x.size());
    }
}

#endif // LBM_GEOM_JET