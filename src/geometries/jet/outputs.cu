#include "outputs.cuh"

#include <vector>
#include <fstream>
#include <cstdint>

#include "properties.cuh"
#include "../../lbm/stencil_active.cuh"

namespace JET
{
    static std::string jet_sections_path(const std::string &out_dir)
    {
        return out_dir + "/outputs/jet_sections.bin";
    }

    static std::string jet_centerline_path(const std::string &out_dir)
    {
        return out_dir + "/outputs/jet_centerline.bin";
    }

    static void write_jet_sections(const LBMState &S, int t, const std::string &out_dir, const DomainTags &tags)
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

                const real_t ux = (host_field<MomentId::ux>(S)[idx] / Stencil::as2 + host_field<MomentId::ux>(S)[idx_n_safe] / Stencil::as2) * real_t(0.5);
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

    static void write_jet_centerline(const LBMState &S, int t, const std::string &out_dir, const DomainTags &tags)
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

            const real_t ux = (host_field<MomentId::ux>(S)[idx] / Stencil::as2 + host_field<MomentId::ux>(S)[idx_n_safe] / Stencil::as2) * real_t(0.5);
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

    void outputs(const LBMState &state, int t, const std::string &out_dir, const DomainTags &tags)
    {
        write_jet_sections(state, t, out_dir, tags);
        write_jet_centerline(state, t, out_dir, tags);
    }
}