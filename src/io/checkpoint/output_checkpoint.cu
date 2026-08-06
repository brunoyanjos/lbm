#include "output_checkpoint.cuh"

#include "../../core/cuda_utils.cuh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sstream>

namespace io
{
    namespace
    {
        void write_field(std::ofstream &file, const real_t *field, std::size_t bytes)
        {
            file.write(reinterpret_cast<const char *>(field), static_cast<std::streamsize>(bytes));
        }

        void read_field(std::ifstream &file, real_t *field, std::size_t bytes)
        {
            file.read(reinterpret_cast<char *>(field), static_cast<std::streamsize>(bytes));
        }

        std::filesystem::path resolve_checkpoint_path(const std::string &checkpoint_path_or_dir)
        {
            namespace fs = std::filesystem;

            const fs::path path(checkpoint_path_or_dir);
            if (fs::is_regular_file(path))
                return path;

            if (!fs::is_directory(path))
                throw std::runtime_error("Checkpoint path does not exist: " + path.string());

            fs::path latest;
            for (const auto &entry : fs::directory_iterator(path))
            {
                if (!entry.is_regular_file())
                    continue;

                const fs::path candidate = entry.path();
                const std::string name = candidate.filename().string();
                if (name.rfind("checkpoint", 0) != 0 || candidate.extension() != ".bin")
                    continue;

                if (latest.empty() || candidate.filename().string() > latest.filename().string())
                    latest = candidate;
            }

            if (latest.empty())
                throw std::runtime_error("No checkpoint*.bin files found in: " + path.string());

            return latest;
        }

        void validate_checkpoint_config(const CheckpointConfig &cfg)
        {
            const CheckpointConfig expected = make_checkpoint_config(cfg.step, cfg.cur);

            if (std::memcmp(cfg.magic, expected.magic, sizeof(cfg.magic)) != 0)
                throw std::runtime_error("Invalid checkpoint magic");
            if (cfg.version != expected.version)
                throw std::runtime_error("Unsupported checkpoint version");
            if (cfg.header_bytes != sizeof(CheckpointConfig))
                throw std::runtime_error("Checkpoint header size does not match this executable");
            if (cfg.endian_marker != expected.endian_marker)
                throw std::runtime_error("Checkpoint endian marker mismatch");
            if (cfg.nx != expected.nx || cfg.ny != expected.ny)
                throw std::runtime_error("Checkpoint grid does not match this executable");
            if (cfg.q != expected.q || cfg.stencil != expected.stencil)
                throw std::runtime_error("Checkpoint stencil does not match this executable");
            if (cfg.real_bytes != expected.real_bytes || cfg.real_is_double != expected.real_is_double)
                throw std::runtime_error("Checkpoint real_t precision does not match this executable");
            if (cfg.reg_order != expected.reg_order || cfg.recurrence != expected.recurrence)
                throw std::runtime_error("Checkpoint regularization config does not match this executable");
            if (cfg.field_count != expected.field_count ||
                cfg.node_count != expected.node_count ||
                cfg.field_bytes != expected.field_bytes ||
                cfg.payload_bytes != expected.payload_bytes)
                throw std::runtime_error("Checkpoint payload layout does not match this executable");
            if (cfg.cur != 0 && cfg.cur != 1)
                throw std::runtime_error("Checkpoint cur buffer must be 0 or 1");
        }

        template <int RegOrder, bool Rec, bool HighOrder>
        void write_extra_fields(std::ofstream &, const LBMStateFor<RegOrder, Rec, HighOrder> &)
        {
        }

        template <bool HighOrder>
        void write_extra_fields(std::ofstream &file, const LBMStateFor<3, false, HighOrder> &state)
        {
            if constexpr (HighOrder)
                write_field(file, state.h_mxxx, state.bytes_field);
            write_field(file, state.h_mxxy, state.bytes_field);
            write_field(file, state.h_mxyy, state.bytes_field);
            if constexpr (HighOrder)
                write_field(file, state.h_myyy, state.bytes_field);
        }

        template <int RegOrder, bool Rec, bool HighOrder>
        void read_extra_fields(std::ifstream &, LBMStateFor<RegOrder, Rec, HighOrder> &)
        {
        }

        template <bool HighOrder>
        void read_extra_fields(std::ifstream &file, LBMStateFor<3, false, HighOrder> &state)
        {
            if constexpr (HighOrder)
                read_field(file, state.h_mxxx, state.bytes_field);
            read_field(file, state.h_mxxy, state.bytes_field);
            read_field(file, state.h_mxyy, state.bytes_field);
            if constexpr (HighOrder)
                read_field(file, state.h_myyy, state.bytes_field);
        }

        template <int RegOrder, bool Rec, bool HighOrder>
        void upload_extra_fields(const LBMStateFor<RegOrder, Rec, HighOrder> &, int)
        {
        }

        template <bool HighOrder>
        void upload_extra_fields(const LBMStateFor<3, false, HighOrder> &state, int c)
        {
            if constexpr (HighOrder)
                CUDA_CHECK(cudaMemcpy(state.d_mxxx[c], state.h_mxxx, state.bytes_field, cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(state.d_mxxy[c], state.h_mxxy, state.bytes_field, cudaMemcpyHostToDevice));
            CUDA_CHECK(cudaMemcpy(state.d_mxyy[c], state.h_mxyy, state.bytes_field, cudaMemcpyHostToDevice));
            if constexpr (HighOrder)
                CUDA_CHECK(cudaMemcpy(state.d_myyy[c], state.h_myyy, state.bytes_field, cudaMemcpyHostToDevice));
        }
    }

    void write_checkpoint_current(const LBMState &state, int step, const std::string &out_dir)
    {
        namespace fs = std::filesystem;

        fs::path checkpoint_dir = fs::path(out_dir) / "checkpoints";
        fs::create_directories(checkpoint_dir);

        const int t_star = step / SAVE_INTERVAL;
        std::ostringstream filename;
        filename << "checkpoint_tstar_" << std::setw(7) << std::setfill('0') << t_star
                 << "_step_" << std::setw(9) << std::setfill('0') << step << ".bin";

        const fs::path filepath = checkpoint_dir / filename.str();
        std::ofstream file(filepath, std::ios::binary);
        if (!file.is_open())
        {
            std::cerr << "Could not open checkpoint file for writing: " << filepath.string() << "\n";
            return;
        }

        const CheckpointConfig cfg = make_checkpoint_config(step, state.cur);
        file.write(reinterpret_cast<const char *>(&cfg), sizeof(cfg));

        write_field(file, state.h_rhoA, state.bytes_field);
        write_field(file, state.h_uxA, state.bytes_field);
        write_field(file, state.h_uyA, state.bytes_field);
        write_field(file, state.h_mxxA, state.bytes_field);
        write_field(file, state.h_mxyA, state.bytes_field);
        write_field(file, state.h_myyA, state.bytes_field);
        write_extra_fields(file, state);

        if (!file.good())
            std::cerr << "Checkpoint write failed: " << filepath.string() << "\n";
    }

    CheckpointConfig read_checkpoint_current(LBMState &state, const std::string &checkpoint_path_or_dir)
    {
        const std::filesystem::path filepath = resolve_checkpoint_path(checkpoint_path_or_dir);
        std::ifstream file(filepath, std::ios::binary);
        if (!file.is_open())
            throw std::runtime_error("Could not open checkpoint file for reading: " + filepath.string());

        CheckpointConfig cfg{};
        file.read(reinterpret_cast<char *>(&cfg), sizeof(cfg));
        if (!file.good())
            throw std::runtime_error("Could not read checkpoint header: " + filepath.string());

        validate_checkpoint_config(cfg);

        state.cur = cfg.cur;
        read_field(file, state.h_rhoA, state.bytes_field);
        read_field(file, state.h_uxA, state.bytes_field);
        read_field(file, state.h_uyA, state.bytes_field);
        read_field(file, state.h_mxxA, state.bytes_field);
        read_field(file, state.h_mxyA, state.bytes_field);
        read_field(file, state.h_myyA, state.bytes_field);
        read_extra_fields(file, state);

        if (!file.good())
            throw std::runtime_error("Could not read checkpoint payload: " + filepath.string());

        const int c = state.cur;
        CUDA_CHECK(cudaMemcpy(state.d_rhoA[c], state.h_rhoA, state.bytes_field, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(state.d_uxA[c], state.h_uxA, state.bytes_field, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(state.d_uyA[c], state.h_uyA, state.bytes_field, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(state.d_mxxA[c], state.h_mxxA, state.bytes_field, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(state.d_mxyA[c], state.h_mxyA, state.bytes_field, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(state.d_myyA[c], state.h_myyA, state.bytes_field, cudaMemcpyHostToDevice));
        upload_extra_fields(state, c);

        std::cout << "[CHECKPOINT] loaded " << filepath.string()
                  << " at step " << cfg.step << "\n";

        return cfg;
    }
} // namespace io
