#include "simulation.cuh"

#include "benchmark.cuh"
#include "progress.cuh"

#include "../io/output_meta.cuh"
#include "../io/output_vtk.cuh"
#include "../io/output_tags_vtk.cuh"

#include "../lbm/state/lbm_state.cuh"
#include "../lbm/lbm_init_state.cuh"
#include "../lbm/lbm_mom_step.cuh"
#include "../lbm/domain/cavity_square_tags.cuh"
#include "../lbm/domain/tags_debug.cuh"

#include "../core/cuda_utils.cuh"
#include "../core/simulation_config.h"
#include "../core/geometry.h"

#include <chrono>
#include <iostream>

namespace app
{
    void run(const CudaConfig &cfg, const RunContext &ctx)
    {
        auto state = lbm_allocate_state();
        init_state(state, cfg);

        DomainTags tags = domain_tags_allocate(ctx.verbose);
        build_cavity_square_tags(tags);

        if (ctx.verbose)
        {
            validate_cavity_square_tags_host(tags, true);

            // evita misturar com a barra (stderr)
            if (ctx.show_progress)
                progress::ProgressUI::suspend_for_log();
            io::write_tags_vtk(tags, ctx.out_dir);
        }

        // ---------------- warmup ----------------
        for (int t = 0; t < ctx.warmup_steps; ++t)
        {
            lbm_mom_step(state, cfg, tags);
            state.cur ^= 1;
        }
        CUDA_CHECK(cudaDeviceSynchronize());

        // ---------------- timers ----------------
        using clock = std::chrono::steady_clock;
        const auto wall0 = clock::now();

        const int t_begin = ctx.warmup_steps;
        const int t_end = N_STEPS;

        // progresso (UI)
        progress::ProgressUI ui;
        ui.start(t_begin, t_end, ctx.show_progress, ctx.progress_hz);

        // timer GPU total
        GpuTimer gt;
        gt.start();

        // eventos para medir GPU parcial sem custo alto
        cudaEvent_t ev_prog0 = nullptr, ev_prog1 = nullptr;
        if (ctx.show_progress)
        {
            CUDA_CHECK(cudaEventCreate(&ev_prog0));
            CUDA_CHECK(cudaEventCreate(&ev_prog1));
            CUDA_CHECK(cudaEventRecord(ev_prog0));
        }

        // ---------------- main loop (mede) ----------------
        for (int t = t_begin; t < t_end; ++t)
        {
            lbm_mom_step(state, cfg, tags);
            state.cur ^= 1;

            // IO (pode ter prints internos) -> limpa a barra antes
            if (ctx.enable_io && (t % SAVE_INTERVAL == 0))
            {
                upload_state_to_host(state);
                io::write_vtk(state, cfg, t, ctx.out_dir);
            }

            // barra (limite de frequência dentro do ProgressUI)
            if (ctx.show_progress && (ui.should_print() || t == t_end - 1))
            {
                const auto now = clock::now();

                // mede GPU parcial desde ev_prog0
                CUDA_CHECK(cudaEventRecord(ev_prog1));
                CUDA_CHECK(cudaEventSynchronize(ev_prog1));

                float ms = 0.0f;
                CUDA_CHECK(cudaEventElapsedTime(&ms, ev_prog0, ev_prog1));
                const double gpu_elapsed_s = double(ms) * 1e-3;

                const double wall_elapsed_s = std::chrono::duration<double>(now - wall0).count();

                const int done_steps = (t - t_begin + 1);
                const double updates = double(NX) * double(NY) * double(done_steps);
                const double mlups_partial = (gpu_elapsed_s > 0.0) ? (updates / gpu_elapsed_s / 1e6) : 0.0;

                ui.print(t, wall_elapsed_s, gpu_elapsed_s, mlups_partial);
            }
        }

        const double gpu_s = gt.stop_seconds();
        const auto wall1 = clock::now();
        const double wall_s = std::chrono::duration<double>(wall1 - wall0).count();

        // ---------------- metrics ----------------
        BenchmarkResult r;
        r.gpu_seconds = gpu_s;
        r.wall_seconds = wall_s;
        r.measured_steps = (N_STEPS - ctx.warmup_steps);

        const double updates = double(NX) * double(NY) * double(r.measured_steps);
        r.mlups_gpu = updates / r.gpu_seconds / 1e6;
        r.mlups_wall = updates / r.wall_seconds / 1e6;

        io::write_performance(ctx.out_dir, cfg, r);

        ui.finish(false);

        if (ctx.show_progress)
        {
            CUDA_CHECK(cudaEventDestroy(ev_prog0));
            CUDA_CHECK(cudaEventDestroy(ev_prog1));
        }

        std::cout << "Simulation finished after " << N_STEPS << " timesteps.\n";
        std::cout << "GPU(s)=" << r.gpu_seconds << " WALL(s)=" << r.wall_seconds
                  << " MLUPS_GPU=" << r.mlups_gpu << " MLUPS_WALL=" << r.mlups_wall << "\n";

        domain_tags_free(tags);
        lbm_free_state(state);
    }
} // namespace app
