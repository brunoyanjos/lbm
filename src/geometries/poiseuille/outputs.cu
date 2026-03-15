#include "outputs.cuh"

#include <vector>
#include <fstream>
#include <cstdint>

#include "properties.cuh"
#include "../../lbm/stencil_active.cuh"

namespace POISEUILLE
{
    static std::string channel_profile_path(const std::string &out_dir)
    {
        return out_dir + "/outputs/channel_profile.bin";
    }

    static void write_channel_profile(const LBMState &S, int t, const std::string &out_dir, const DomainTags &tags)
    {
        // seção “desenvolvida”
        const int x_sample = NX / 2;

        std::vector<float> ux_y;
        ux_y.resize(NY);

        for (int y = 0; y < NY; ++y)
        {
            const int idx = x_sample + NX * y;
            const int idx_n = (x_sample - 1) + NX * y;

            if (tags.h_node && tags.h_node[idx] == to_u8(NodeId::SOLID))
            {
                ux_y[y] = 0.0f;
                continue;
            }

            const real_t ux = (host_field<MomentId::ux>(S)[idx] / Stencil::as2 + host_field<MomentId::ux>(S)[idx_n] / Stencil::as2) * real_t(0.5);
            ux_y[y] = (float)ux;
        }

        const std::string path = channel_profile_path(out_dir);
        std::ofstream f(path, std::ios::binary);

        const char magic[4] = {'C', 'H', 'P', '1'};
        f.write(magic, 4);

        const int32_t nx32 = NX;
        const int32_t ny32 = NY;
        const int32_t t32 = (int32_t)t;
        const int32_t xs32 = (int32_t)x_sample;

        f.write(reinterpret_cast<const char *>(&nx32), sizeof(nx32));
        f.write(reinterpret_cast<const char *>(&ny32), sizeof(ny32));
        f.write(reinterpret_cast<const char *>(&t32), sizeof(t32));
        f.write(reinterpret_cast<const char *>(&xs32), sizeof(xs32));

        f.write(reinterpret_cast<const char *>(ux_y.data()), sizeof(float) * ux_y.size());
    }

    void outputs(const LBMState &state, int t, const std::string &out_dir, const DomainTags &tags)
    {
        write_channel_profile(state, t, out_dir, tags);
    }
}