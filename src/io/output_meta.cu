#include "output_meta.cuh"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include "../core/physics.h"
#include "../core/geometry.h"

namespace io
{
    void write_performance(const std::string &out_dir,
                           const CudaConfig &,
                           const app::BenchmarkResult &r)
    {
        namespace fs = std::filesystem;
        fs::create_directories(fs::path(out_dir) / "meta");

        std::ofstream f(fs::path(out_dir) / "meta" / "performance.txt");
        f << std::fixed << std::setprecision(6);
        f << "NX " << NX << "\n";
        f << "NY " << NY << "\n";
        f << "measured_steps " << r.measured_steps << "\n";
        f << "gpu_seconds " << r.gpu_seconds << "\n";
        f << "wall_seconds " << r.wall_seconds << "\n";
        f << "MLUPS_GPU " << r.mlups_gpu << "\n";
        f << "MLUPS_Wall " << r.mlups_wall << "\n";
    }
}
