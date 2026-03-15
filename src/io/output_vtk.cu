#include "output_vtk.cuh"

#include "../geometries/active_geometry.cuh"
#include "../core/indexing.cuh"
#include "../lbm/stencil_active.cuh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

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

        const int nx = Geometry::NX;
        const int ny = Geometry::NY;

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
                const real_t rho = host_field<MomentId::rho>(S)[idx] + Geometry::RHO_0;
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
                const real_t ux = host_field<MomentId::ux>(S)[idx] / Stencil::as2;
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
                const real_t uy = host_field<MomentId::uy>(S)[idx] / Stencil::as2;
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
                const real_t ux = host_field<MomentId::ux>(S)[idx] / Stencil::as2;
                const real_t uy = host_field<MomentId::uy>(S)[idx] / Stencil::as2;

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

} // namespace io
