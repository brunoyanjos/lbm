#include "output_vtk.cuh"

#include "../../core/geometry.h"
#include "../../core/physics.h"
#include "../../core/indexing.cuh"
#include "../../lbm/stencil_active.cuh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>

namespace io
{
    __host__ void write_vti(const LBMState &S, const CudaConfig & /*cfg*/, int step, const std::string &out_dir)
    {
        namespace fs = std::filesystem;

        fs::path vtk_dir = fs::path(out_dir) / "vtk";
        fs::create_directories(vtk_dir);

        std::ostringstream filename;
        filename << "output_" << std::setw(6) << std::setfill('0') << step << ".vti";

        fs::path filepath = vtk_dir / filename.str();

        std::ofstream file(filepath.string());
        if (!file.is_open())
        {
            std::cerr << "Could not open VTI file for writing: " << filepath.string() << "\n";
            return;
        }

        const int nx = NX;
        const int ny = NY;

        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        file << "  <ImageData WholeExtent=\"0 " << (nx - 1)
             << " 0 " << (ny - 1)
             << " 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
        file << "    <Piece Extent=\"0 " << (nx - 1)
             << " 0 " << (ny - 1)
             << " 0 0\">\n";

        file << "      <PointData>\n";

        // rho
        file << "        <DataArray type=\"Float32\" Name=\"rho\" format=\"ascii\">\n";
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const size_t idx = idxGlobal(x, y);
                const real_t rho = S.h_rho[idx] + RHO_0;
                file << "          " << static_cast<float>(rho) << "\n";
            }
        }
        file << "        </DataArray>\n";

        // ux
        file << "        <DataArray type=\"Float32\" Name=\"ux\" format=\"ascii\">\n";
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const size_t idx = idxGlobal(x, y);
                const real_t ux = S.h_ux[idx] / Stencil::as2;
                file << "          " << static_cast<float>(ux) << "\n";
            }
        }
        file << "        </DataArray>\n";

        // uy
        file << "        <DataArray type=\"Float32\" Name=\"uy\" format=\"ascii\">\n";
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const size_t idx = idxGlobal(x, y);
                const real_t uy = S.h_uy[idx] / Stencil::as2;
                file << "          " << static_cast<float>(uy) << "\n";
            }
        }
        file << "        </DataArray>\n";

        // velocity
        file << "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const size_t idx = idxGlobal(x, y);
                const real_t ux = S.h_ux[idx] / Stencil::as2;
                const real_t uy = S.h_uy[idx] / Stencil::as2;

                file << "          "
                     << static_cast<float>(ux) << " "
                     << static_cast<float>(uy) << " "
                     << 0.0f << "\n";
            }
        }
        file << "        </DataArray>\n";

        file << "      </PointData>\n";
        file << "      <CellData>\n";
        file << "      </CellData>\n";
        file << "    </Piece>\n";
        file << "  </ImageData>\n";
        file << "</VTKFile>\n";

        file.close();
    }

    namespace
    {
        void copy_local_interior_to_global(real_t *dst, const real_t *src, const LocalDomain &domain)
        {
            for (int y_local = 0; y_local < domain.local_ny; ++y_local)
            {
                const int y_global = domain.y_begin + y_local;
                for (int x = 0; x < domain.nx; ++x)
                {
                    const size_t src_idx = idxLocal(x, y_local + domain.halo, domain.nx);
                    const size_t dst_idx = idxGlobal(x, y_global);
                    dst[dst_idx] = src[src_idx];
                }
            }
        }
    }

    __host__ void write_vti(const std::vector<LBMState> &states, int step, const std::string &out_dir)
    {
        LBMState snapshot = lbm_allocate_state(make_local_domain(0, NY, 0));

        for (const LBMState &state : states)
        {
            copy_local_interior_to_global(snapshot.h_rho, state.h_rho, state.domain);
            copy_local_interior_to_global(snapshot.h_ux, state.h_ux, state.domain);
            copy_local_interior_to_global(snapshot.h_uy, state.h_uy, state.domain);
            copy_local_interior_to_global(snapshot.h_mxx, state.h_mxx, state.domain);
            copy_local_interior_to_global(snapshot.h_mxy, state.h_mxy, state.domain);
            copy_local_interior_to_global(snapshot.h_myy, state.h_myy, state.domain);
        }

        write_vti(snapshot, make_config(NX, NY), step, out_dir);
        lbm_free_state(snapshot);
    }

}
