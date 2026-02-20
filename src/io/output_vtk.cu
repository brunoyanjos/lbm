#include "output_vtk.cuh"

#include "../core/geometry.h"
#include "../core/physics.h"
#include "../core/indexing.cuh"
#include "../lbm/stencil_active.cuh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

namespace io
{
    __host__ void write_vtk(const LBMState &S, const CudaConfig & /*cfg*/, int step, const std::string &out_dir)
    {
        namespace fs = std::filesystem;

        fs::path vtk_dir = fs::path(out_dir) / "vtk";
        fs::create_directories(vtk_dir);

        std::ostringstream filename;
        filename << "output_" << std::setw(6) << std::setfill('0') << step << ".vtk";

        fs::path filepath = vtk_dir / filename.str();

        std::ofstream file(filepath.string());
        if (!file.is_open())
        {
            std::cerr << "Could not open VTK file for writing: " << filepath.string() << "\n";
            return;
        }

        // --- HEADER ---
        file << "# vtk DataFile Version 3.0\n";
        file << "LBM output\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_POINTS\n";
        file << "DIMENSIONS " << NX << " " << NY << " 1\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING 1 1 1\n";
        file << "POINT_DATA " << (static_cast<size_t>(NX) * static_cast<size_t>(NY)) << "\n";

        // --- DENSITY ---
        file << "SCALARS rho float 1\n";
        file << "LOOKUP_TABLE default\n";
        for (int y = 0; y < NY; ++y)
        {
            for (int x = 0; x < NX; ++x)
            {
                const size_t idx = idxGlobal(x, y);
                const real_t rho = S.h_rho[idx] + RHO_0; // h_rho stores (rho - RHO_0)
                file << static_cast<float>(rho) << "\n";
            }
        }

        // --- VELOCITY ---
        file << "VECTORS velocity float\n";
        for (int y = 0; y < NY; ++y)
        {
            for (int x = 0; x < NX; ++x)
            {
                const size_t idx = idxGlobal(x, y);

                // h_ux/h_uy store scaled velocity (as2*u) by current convention
                const real_t ux = S.h_ux[idx] / Stencil::as2;
                const real_t uy = S.h_uy[idx] / Stencil::as2;

                file << static_cast<float>(ux) << " "
                     << static_cast<float>(uy) << " "
                     << 0.0f << "\n";
            }
        }

        file.close();
    }

} // namespace io
