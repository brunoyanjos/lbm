#include "output_annul_profile_bin.cuh"
#include "../../../lbm/stencil_active.cuh"

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include "../../../core/active_geometry.cuh"
#include "../../../core/types.cuh"

namespace io
{
    static std::string annul_profile_path(const std::string &out_dir)
    {
        return out_dir + "/outputs/annul_profile.bin";
    }

    static inline void cart_to_polar(real_t ux, real_t uy, real_t c, real_t s,
                                     real_t &ur, real_t &utheta)
    {
        ur = ux * c + uy * s;
        utheta = -ux * s + uy * c;
    }

    void write_annul_profile(const LBMState &state, int t, const std::string &out_dir, const DomainTags &tags)
    {
        const real_t xc = r_cast(NX - 1) * r::half;
        const real_t yc = r_cast(NY - 1) * r::half;

        const int y_line = NY / 2;

        // varre metade direita do domínio
        const int x_begin = NX / 2;
        const int x_end = NX; // exclusivo

        std::vector<float> r_prof;
        std::vector<float> ur_prof;
        std::vector<float> ut_prof;

        // (opcional) reserva um chute pra evitar realloc
        r_prof.reserve(x_end - x_begin);
        ur_prof.reserve(x_end - x_begin);
        ut_prof.reserve(x_end - x_begin);

        for (int x = x_begin; x < x_end; ++x)
        {
            const int idx = x + NX * y_line;
            const int idx_n = x + NX * (y_line - 1);

            if (tags.h_node[idx] == to_u8(NodeId::SOLID))
                continue;

            if (tags.h_node[idx_n] == to_u8(NodeId::SOLID))
                continue;

            const real_t ux = (state.h_ux[idx] / Stencil::as2 + state.h_ux[idx_n] / Stencil::as2) * r::half;
            const real_t uy = (state.h_uy[idx] / Stencil::as2 + state.h_uy[idx_n] / Stencil::as2) * r::half;

            const real_t dx = r_cast(x) - xc;
            const real_t y_mid = (r_cast(y_line) + r_cast(y_line - 1)) * r::half;
            const real_t dy = y_mid - yc;

            const real_t r2 = dx * dx + dy * dy;
            const real_t invr = rsqrt(r2 + real_t(1e-30));
            const real_t c = dx * invr;
            const real_t s = dy * invr;
            const real_t r = r2 * invr;

            real_t ur, ut;
            cart_to_polar(ux, uy, c, s, ur, ut);

            r_prof.push_back((float)r);
            ur_prof.push_back((float)ur);
            ut_prof.push_back((float)ut);
        }

        const int32_t n32 = (int32_t)r_prof.size();
        const int32_t t32 = (int32_t)t;
        const int32_t y32 = (int32_t)y_line;
        const int32_t x032 = (int32_t)x_begin; // agora é o começo do sweep, não o primeiro fluido

        const std::string path = annul_profile_path(out_dir);
        std::ofstream f(path, std::ios::binary);

        const char magic[4] = {'A', 'P', 'R', '1'};
        f.write(magic, 4);

        f.write(reinterpret_cast<const char *>(&n32), sizeof(n32));
        f.write(reinterpret_cast<const char *>(&t32), sizeof(t32));
        f.write(reinterpret_cast<const char *>(&y32), sizeof(y32));
        f.write(reinterpret_cast<const char *>(&x032), sizeof(x032));

        if (n32 > 0)
        {
            f.write(reinterpret_cast<const char *>(r_prof.data()), sizeof(float) * r_prof.size());
            f.write(reinterpret_cast<const char *>(ur_prof.data()), sizeof(float) * ur_prof.size());
            f.write(reinterpret_cast<const char *>(ut_prof.data()), sizeof(float) * ut_prof.size());
        }
    }
} // namespace io