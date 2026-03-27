#include "debug_domain.cuh"

#include <cstdio>
#include <cstdint>

#include "core/geometry.h"
#include "lbm/domain/build_tags.cuh"

namespace io
{
    void debug_domain(const DomainTags &T,
                      bool print_domain,
                      bool print_boundary_ids,
                      int max_boundary_points)
    {
        if (!T.h_node)
        {
            std::fprintf(stderr,
                         "[debug_domain] T.h_node is null. "
                         "Make sure build_tags() copied back.\n");
            return;
        }

        const uint8_t *nodes = T.h_node;

        const uint8_t FLUID = to_u8(NodeId::FLUID);
        const uint8_t SOLID = to_u8(NodeId::SOLID);

        if (print_domain)
        {
            std::printf("\n=== DOMAIN (NY=%d, NX=%d) ===\n", int(NY), int(NX));
            std::printf("Legend: F=FLUID, S=SOLID, 1..14=BOUNDARY ID\n\n");

            for (int y = int(NY) - 1; y >= 0; --y)
            {
                std::printf("%4d |", y);

                for (int x = 0; x < int(NX); ++x)
                {
                    const size_t idx = size_t(x) + size_t(NX) * size_t(y);
                    const uint8_t nid = nodes[idx];

                    if (nid == SOLID)
                        std::printf("%3s", "S");
                    else if (nid == FLUID)
                        std::printf("%3s", "F");
                    else
                        std::printf("%3u", unsigned(nid));
                }

                std::printf("\n");
            }

            std::printf("      +");
            for (int x = 0; x < int(NX) * 3; ++x)
                std::printf("-");
            std::printf("\n");

            std::printf("       ");
            for (int x = 0; x < int(NX); ++x)
                std::printf("%3d", x % 100);
            std::printf("\n\n");
        }

        if (print_boundary_ids)
        {
            std::printf("=== BOUNDARY IDS ===\n");
            std::printf("Each boundary node prints node_id and bit pattern.\n\n");

            int printed = 0;

            for (int y = 0; y < int(NY); ++y)
            {
                for (int x = 0; x < int(NX); ++x)
                {
                    const size_t idx = size_t(x) + size_t(NX) * size_t(y);
                    const uint8_t nid = nodes[idx];

                    if (nid == SOLID || nid == FLUID)
                        continue;

                    const unsigned b1 = (nid & 1u) ? 1u : 0u;
                    const unsigned b2 = (nid & 2u) ? 1u : 0u;
                    const unsigned b4 = (nid & 4u) ? 1u : 0u;
                    const unsigned b8 = (nid & 8u) ? 1u : 0u;

                    std::printf("(x=%d, y=%d) idx=%zu  node_id=%u  bits=[%u %u %u %u]\n",
                                x, y, idx, unsigned(nid), b1, b2, b4, b8);

                    if (++printed >= max_boundary_points)
                    {
                        std::printf("... truncated (max_boundary_points=%d)\n",
                                    max_boundary_points);
                        std::printf("\n");
                        std::fflush(stdout);
                        return;
                    }
                }
            }

            if (printed == 0)
                std::printf("(no boundary points found)\n");

            std::printf("\n");
        }

        std::fflush(stdout);
    }
}