#include "output_checkpoint.cuh"

#include "checkpoint_config.cuh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace io
{
    namespace
    {
        void write_field(std::ofstream &file, const real_t *field, std::size_t bytes)
        {
            file.write(reinterpret_cast<const char *>(field), static_cast<std::streamsize>(bytes));
        }
    }

    void write_checkpoint_current(const LBMState &state, int step, const std::string &out_dir)
    {
        namespace fs = std::filesystem;

        fs::path checkpoint_dir = fs::path(out_dir) / "checkpoints";
        fs::create_directories(checkpoint_dir);

        std::ostringstream filename;
        filename << "checkpoint_" << std::setw(9) << std::setfill('0') << step << ".bin";

        const fs::path filepath = checkpoint_dir / filename.str();
        std::ofstream file(filepath, std::ios::binary);
        if (!file.is_open())
        {
            std::cerr << "Could not open checkpoint file for writing: " << filepath.string() << "\n";
            return;
        }

        const CheckpointConfig cfg = make_checkpoint_config(step, state.cur);
        file.write(reinterpret_cast<const char *>(&cfg), sizeof(cfg));

        write_field(file, state.h_rho, state.bytes_field);
        write_field(file, state.h_ux, state.bytes_field);
        write_field(file, state.h_uy, state.bytes_field);
        write_field(file, state.h_mxx, state.bytes_field);
        write_field(file, state.h_mxy, state.bytes_field);
        write_field(file, state.h_myy, state.bytes_field);

        if (!file.good())
            std::cerr << "Checkpoint write failed: " << filepath.string() << "\n";
    }
} // namespace io
