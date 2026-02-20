#pragma once
#include "cavity_square_tags.cuh"

#include "../../core/geometry.h"
#include "../../core/indexing.cuh"
#include "../stencil_active.cuh"

#include <cstdint>
#include <iostream>
#include <stdexcept>

static inline int popcount_u32(uint32_t x) { return __builtin_popcount(x); }

static inline uint32_t expected_valid_mask_at(int x, int y)
{
    uint32_t m = 0u;
    m |= (1u << 0); // rest sempre válido

    for (int i = 1; i < Stencil::Q; ++i)
    {
        const int xn = x + Stencil::cx(i);
        const int yn = y + Stencil::cy(i);

        if (xn < 0 || xn >= NX || yn < 0 || yn >= NY)
            continue;
        m |= (1u << i);
    }
    return m;
}

static inline void validate_cavity_square_tags_host(const DomainTags &T, bool strict_masks = true)
{
    if (!T.h_valid || !T.h_wall)
        throw std::runtime_error("validate_cavity_square_tags_host requires host buffers (allocate with host_buffers=true).");

    size_t c_none = 0, c_left = 0, c_right = 0, c_bottom = 0, c_top = 0, c_corner = 0;

    for (int y = 0; y < NY; ++y)
    {
        for (int x = 0; x < NX; ++x)
        {
            const size_t idx = idxGlobal(x, y);
            const uint8_t w = T.h_wall[idx];

            switch ((WallId)w)
            {
            case WallId::NONE:
                c_none++;
                break;
            case WallId::LEFT:
                c_left++;
                break;
            case WallId::RIGHT:
                c_right++;
                break;
            case WallId::BOTTOM:
                c_bottom++;
                break;
            case WallId::TOP:
                c_top++;
                break;
            case WallId::CORNER:
                c_corner++;
                break;
            default:
                throw std::runtime_error("Unknown WallId value found in tag buffer.");
            }
        }
    }

    const size_t exp_corner = 4;
    const size_t exp_left = (NY >= 2) ? size_t(NY - 2) : 0;
    const size_t exp_right = (NY >= 2) ? size_t(NY - 2) : 0;
    const size_t exp_bottom = (NX >= 2) ? size_t(NX - 2) : 0;
    const size_t exp_top = (NX >= 2) ? size_t(NX - 2) : 0;
    const size_t exp_none = (NX >= 2 && NY >= 2) ? size_t(NX - 2) * size_t(NY - 2) : 0;

    std::cout << "[TAGS] counts: NONE=" << c_none
              << " LEFT=" << c_left
              << " RIGHT=" << c_right
              << " BOTTOM=" << c_bottom
              << " TOP=" << c_top
              << " CORNER=" << c_corner
              << " (NX=" << NX << " NY=" << NY << ")\n";

    if (c_corner != exp_corner || c_left != exp_left || c_right != exp_right ||
        c_bottom != exp_bottom || c_top != exp_top || c_none != exp_none)
    {
        std::cerr << "[TAGS] expected: NONE=" << exp_none
                  << " LEFT=" << exp_left
                  << " RIGHT=" << exp_right
                  << " BOTTOM=" << exp_bottom
                  << " TOP=" << exp_top
                  << " CORNER=" << exp_corner << "\n";
        throw std::runtime_error("WallId counts do not match expected cavity-square layout.");
    }

    if (strict_masks)
    {
        int mismatches = 0;
        const int max_print = 16;

        for (int y = 0; y < NY; ++y)
        {
            for (int x = 0; x < NX; ++x)
            {
                const size_t idx = idxGlobal(x, y);
                const uint32_t got = T.h_valid[idx];
                const uint32_t exp = expected_valid_mask_at(x, y);

                if (got != exp)
                {
                    if (mismatches < max_print)
                    {
                        std::cerr << "[TAGS] valid mismatch at (" << x << "," << y << ") "
                                  << "got=0x" << std::hex << got
                                  << " exp=0x" << exp << std::dec
                                  << " got_bits=" << popcount_u32(got)
                                  << " exp_bits=" << popcount_u32(exp)
                                  << "\n";
                    }
                    ++mismatches;
                }
            }
        }

        if (mismatches)
        {
            std::cerr << "[TAGS] total valid-mask mismatches: " << mismatches << "\n";
            throw std::runtime_error("Valid-mask buffer differs from expected stencil-neighborhood mask.");
        }
    }

    // prints úteis (cantos e meio das bordas)
    auto dump = [&](int x, int y)
    {
        const size_t idx = idxGlobal(x, y);
        const uint32_t m = T.h_valid[idx];
        const uint8_t w = T.h_wall[idx];
        std::cout << "[TAGS] (" << x << "," << y << ") wall=" << int(w)
                  << " valid=0x" << std::hex << m << std::dec
                  << " bits=" << popcount_u32(m) << "\n";
    };

    dump(0, 0);
    dump(NX - 1, 0);
    dump(0, NY - 1);
    dump(NX - 1, NY - 1);
    dump(NX / 2, 0);
    dump(NX / 2, NY - 1);
    dump(0, NY / 2);
    dump(NX - 1, NY / 2);

    std::cout << "[TAGS] validation OK.\n";
}
