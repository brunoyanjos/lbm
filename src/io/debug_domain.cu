#include "debug_domain.cuh"

#include <cstdio>
#include <cstdint>
#include <algorithm>

#include "lbm/stencil_active.cuh" // Stencil::Q, cx, cy
#include "core/types.cuh"         // real_t, etc (se precisar)
#include "core/cuda_utils.cuh"    // CUDA_CHECK (se quiser)
#include "core/math_utils.cuh"    // CUDA_CHECK (se quiser)
#include "core/geometry.h"

namespace
{
    // full mask para Q <= 32
    static inline constexpr mask_t full_mask_host()
    {
        static_assert(Stencil::Q <= 64, "Stencil::Q must be <= 64 for mask_t");

        if constexpr (Stencil::Q == int(sizeof(mask_t) * 8))
            return ~mask_t(0);
        else
            return (mask_t(1) << Stencil::Q) - mask_t(1);
    }

    static inline bool in_bounds(int x, int y)
    {
        return (0 <= x && x < int(NX) && 0 <= y && y < int(NY));
    }

    // acesso host para o node buffer (uint8), fora -> SOLID (conservador)
    static inline uint8_t get_node_safe_host(const uint8_t *nodes, int x, int y, uint8_t SOLID)
    {
        if (!in_bounds(x, y))
            return SOLID;
        return nodes[size_t(x) + size_t(NX) * size_t(y)];
    }

    // detecta se é FLUID perto de SOLID (qualquer vizinho sólido em 1..Q-1)
    static inline bool fluid_near_solid(const uint8_t *nodes, int x, int y, uint8_t SOLID)
    {
        const uint8_t me = get_node_safe_host(nodes, x, y, SOLID);

        for (int i = 1; i < Stencil::Q; ++i)
        {
            const int xn = x + Stencil::cx(i);
            const int yn = y + Stencil::cy(i);
            if (get_node_safe_host(nodes, xn, yn, SOLID) == SOLID)
                return true;
        }
        return false;
    }

    static inline char node_char(uint8_t node_id, uint8_t FLUID, uint8_t SOLID, uint8_t DIRICHLET, uint8_t INLET, uint8_t OUTLET)
    {
        if (node_id == SOLID)
            return 'S';
        if (node_id == DIRICHLET)
            return 'D';
        if (node_id == INLET)
            return 'I';
        if (node_id == OUTLET)
            return 'O';
        if (node_id == FLUID)
            return 'F';
        return '?';
    }
}

namespace io
{
    void debug_domain(const DomainTags &T,
                      bool print_domain,
                      bool print_masks,
                      int max_mask_points)
    {
        if (!T.h_node || !T.h_valid)
        {
            std::fprintf(stderr,
                         "[debug_domain] T.h_node/T.h_valid are null. "
                         "Enable host buffers and make sure build_tags() copied back.\n");
            return;
        }

        const uint8_t *nodes = T.h_node;
        const mask_t *valid = T.h_valid;

        // ajuste esses casts/ids conforme seu enum NodeId
        const uint8_t FLUID = to_u8(NodeId::FLUID);
        const uint8_t SOLID = to_u8(NodeId::SOLID);
        const uint8_t DIRICHLET = to_u8(NodeId::DIRICHLET);
        const uint8_t INLET = to_u8(NodeId::INLET);
        const uint8_t OUTLET = to_u8(NodeId::OUTLET);

        const mask_t FM = full_mask_host();

        if (print_domain)
        {
            std::printf("\n=== DOMAIN (NY=%d, NX=%d) ===\n", int(NY), int(NX));
            std::printf("Legend: F=FLUID, S=SOLID, D=DIRICHLET\n\n");

            // imprime com y decrescente para o "topo" aparecer em cima
            for (int y = int(NY) - 1; y >= 0; --y)
            {
                std::printf("%4d | ", y);
                for (int x = 0; x < int(NX); ++x)
                {
                    const size_t idx = size_t(x) + size_t(NX) * size_t(y);
                    const char c = node_char(nodes[idx], FLUID, SOLID, DIRICHLET, INLET, OUTLET);
                    std::printf("%c", c);
                }
                std::printf("\n");
            }

            std::printf("      +");
            for (int x = 0; x < int(NX) + 1; ++x)
                std::printf("-");
            std::printf("\n       ");
            for (int x = 0; x < int(NX); ++x)
                std::printf("%d", (x % 10));
            std::printf("\n\n");
        }

        if (print_masks)
        {
            std::printf("=== PARTIAL VALID MASKS (fluid near solid, excluding dirichlet) ===\n");
            std::printf("Mask bits order: i=0..Q-1 (Q=%d)\n\n", int(Stencil::Q));

            int printed = 0;

            // varre tudo, mas imprime só casos relevantes
            for (int y = 0; y < int(NY); ++y)
            {
                for (int x = 0; x < int(NX); ++x)
                {
                    const size_t idx = size_t(x) + size_t(NX) * size_t(y);
                    const uint8_t nid = nodes[idx];

                    if (nid == SOLID)
                        continue;

                    const mask_t m = valid[idx];
                    if (m == FM)
                        continue;

                    if (!fluid_near_solid(nodes, x, y, SOLID))
                        continue;

                    std::printf("(x=%d, y=%d) idx=%zu  mask=0x%016llX  bits: ",
                                x, y, idx, (unsigned long long)m);

                    for (int i = 0; i < Stencil::Q; ++i)
                    {
                        const mask_t bit = (m >> i) & mask_t(1);
                        std::printf("%u", unsigned(bit));
                        if (i + 1 < Stencil::Q)
                            std::printf(" ");
                    }
                    std::printf("\n");

                    if (++printed >= std::max(0, max_mask_points))
                    {
                        std::printf("... truncated (max_mask_points=%d)\n", max_mask_points);
                        y = int(NY); // break duplo
                        break;
                    }
                }
            }

            if (printed == 0)
                std::printf("(no points matched criteria)\n");

            std::printf("\n");
        }

        std::fflush(stdout);
    }
}