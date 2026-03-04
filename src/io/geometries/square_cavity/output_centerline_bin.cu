#include "output_centerline_bin.cuh"
#include "../../../lbm/stencil_active.cuh"

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include "../../../core/active_geometry.cuh"
#include "../../../core/types.cuh"

namespace io
{
    static std::string centerline_path(const std::string &out_dir)
    {
        return out_dir + "/outputs/centerline.bin";
    }

    void write_centerline_profiles(const LBMState &state, int t, const std::string &out_dir)
    {
        const int xc = NX / 2;
        const int yc = NY / 2;

        std::vector<real_t> ux_xc_y(NY);
        std::vector<real_t> uy_yc_x(NX);

        for (int y = 0; y < NY; ++y)
        {
            const int idx = xc + NX * y;
            const int idx_n = (xc - 1) + NX * y;

            ux_xc_y[y] = (state.h_ux[idx] / Stencil::as2 + state.h_ux[idx_n] / Stencil::as2) * real_t(0.5);
        }

        for (int x = 0; x < NX; ++x)
        {
            const int idx = x + NX * yc;
            const int idx_n = x + NX * (yc - 1);

            uy_yc_x[x] = (state.h_uy[idx] / Stencil::as2 + state.h_uy[idx_n] / Stencil::as2) * real_t(0.5);
        }

        const std::string path = centerline_path(out_dir);
        std::ofstream f(path, std::ios::binary);

        // header
        const char magic[4] = {'C', 'L', 'N', '1'};
        f.write(magic, 4);

        const int32_t nx32 = NX;
        const int32_t ny32 = NY;
        const int32_t xc32 = xc;
        const int32_t yc32 = yc;
        const int32_t t32 = (int32_t)t;

        f.write(reinterpret_cast<const char *>(&nx32), sizeof(nx32));
        f.write(reinterpret_cast<const char *>(&ny32), sizeof(ny32));
        f.write(reinterpret_cast<const char *>(&xc32), sizeof(xc32));
        f.write(reinterpret_cast<const char *>(&yc32), sizeof(yc32));
        f.write(reinterpret_cast<const char *>(&t32), sizeof(t32));

        // arrays
        f.write(reinterpret_cast<const char *>(ux_xc_y.data()), sizeof(real_t) * ux_xc_y.size());
        f.write(reinterpret_cast<const char *>(uy_yc_x.data()), sizeof(real_t) * uy_yc_x.size());
    }
} // namespace io