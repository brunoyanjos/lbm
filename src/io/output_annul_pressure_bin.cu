#include "output_annul_pressure_bin.cuh"

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "../core/geometry.h"
#include "../core/types.cuh"
#include "../core/physics.h"

namespace io
{
    static std::string annul_pressure_path(const std::string &out_dir)
    {
        return out_dir + "/outputs/annul_pressure.bin";
    }

    void write_annul_pressure(const LBMState &state, int t, const std::string &out_dir, const DomainTags &tags)
    {
        const real_t xc = r_cast(NX - 1) * r::half;
        const real_t yc = r_cast(NY - 1) * r::half;

        const int y_line = NY / 2;

        // novo sweep: metade direita
        const int x_begin = NX / 2;
        const int x_end = NX; // exclusivo

        std::vector<float> r_prof;
        std::vector<float> rho_p_prof;

        r_prof.reserve(x_end - x_begin);
        rho_p_prof.reserve(x_end - x_begin);

        for (int x = x_begin; x < x_end; ++x)
        {
            const int idx = x + NX * y_line;
            const int idx_n = x + NX * (y_line - 1);

            // filtra sólidos (ajuste se seu array/semântica for diferente)
            if (tags.h_node[idx] == to_u8(NodeId::SOLID))
                continue;
            if (tags.h_node[idx_n] == to_u8(NodeId::SOLID))
                continue;

            // rho_p: aqui você está gravando rho "descontado" (rho - RHO_0),
            // mantendo consistente com o teu state.h_rho.
            // Se você quiser o rho absoluto no bin, é aqui que soma RHO_0.
            const real_t rho_p = (state.h_rho[idx] + state.h_rho[idx_n]) * r::half;

            const real_t dx = r_cast(x) - xc;
            const real_t dy = r_cast(yc) - yc;

            const real_t r2 = dx * dx + dy * dy;
            const real_t invr = rsqrt(r2 + real_t(1e-30));
            const real_t r = r2 * invr;

            r_prof.push_back((float)r);
            rho_p_prof.push_back((float)rho_p);
        }

        const std::string path = annul_pressure_path(out_dir);
        std::ofstream f(path, std::ios::binary);

        const char magic[4] = {'A', 'P', 'P', '1'};
        f.write(magic, 4);

        const int32_t n32 = (int32_t)r_prof.size();
        const int32_t t32 = (int32_t)t;
        const int32_t y32 = (int32_t)y_line;
        const int32_t x032 = (int32_t)x_begin;

        f.write(reinterpret_cast<const char *>(&n32), sizeof(n32));
        f.write(reinterpret_cast<const char *>(&t32), sizeof(t32));
        f.write(reinterpret_cast<const char *>(&y32), sizeof(y32));
        f.write(reinterpret_cast<const char *>(&x032), sizeof(x032));

        if (n32 > 0)
        {
            f.write(reinterpret_cast<const char *>(r_prof.data()), sizeof(float) * r_prof.size());
            f.write(reinterpret_cast<const char *>(rho_p_prof.data()), sizeof(float) * rho_p_prof.size());
        }
    }
}