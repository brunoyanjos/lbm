#include "output_tags_vtk.cuh"

#include "../core/geometry.h"
#include "../core/indexing.cuh"
#include "../lbm/stencil_active.cuh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace io
{

    static inline int popcount_u32(uint32_t x) { return __builtin_popcount(x); }

    void write_tags_vtk(const DomainTags &T, const std::string &out_dir)
    {
        if (!T.h_valid || !T.h_wall)
        {
            std::cerr << "[VTK-TAGS] Host buffers missing. Allocate DomainTags with host_buffers=true.\n";
            return;
        }

        namespace fs = std::filesystem;
        fs::path vtk_dir = fs::path(out_dir) / "vtk";
        fs::create_directories(vtk_dir);

        fs::path filepath = vtk_dir / "tags.vtk";

        std::ofstream file(filepath.string());
        if (!file.is_open())
        {
            std::cerr << "[VTK-TAGS] Could not open file: " << filepath.string() << "\n";
            return;
        }

        file << "# vtk DataFile Version 3.0\n";
        file << "LBM tags debug\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_POINTS\n";
        file << "DIMENSIONS " << NX << " " << NY << " 1\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING 1 1 1\n";
        file << "POINT_DATA " << (static_cast<size_t>(NX) * static_cast<size_t>(NY)) << "\n";

        // wall_id
        file << "SCALARS wall_id unsigned_char 1\n";
        file << "LOOKUP_TABLE default\n";
        for (int y = 0; y < NY; ++y)
        {
            for (int x = 0; x < NX; ++x)
            {
                const size_t idx = idxGlobal(x, y);
                file << static_cast<unsigned int>(T.h_wall[idx]) << "\n";
            }
        }

        // valid_bits
        file << "SCALARS valid_bits int 1\n";
        file << "LOOKUP_TABLE default\n";
        for (int y = 0; y < NY; ++y)
        {
            for (int x = 0; x < NX; ++x)
            {
                const size_t idx = idxGlobal(x, y);
                const uint32_t m = T.h_valid[idx];
                file << popcount_u32(m) << "\n";
            }
        }

        // valid_mask (como unsigned int)
        file << "SCALARS valid_mask unsigned_int 1\n";
        file << "LOOKUP_TABLE default\n";
        for (int y = 0; y < NY; ++y)
        {
            for (int x = 0; x < NX; ++x)
            {
                const size_t idx = idxGlobal(x, y);
                file << static_cast<uint32_t>(T.h_valid[idx]) << "\n";
            }
        }

        file.close();
        std::cout << "[VTK-TAGS] wrote: " << filepath.string() << "\n";
    }

} // namespace io
