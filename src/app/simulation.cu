#include "simulation.cuh"

#include "benchmark.cuh"
#include "progress.cuh"

#include "../io/meta/output_meta.cuh"
#include "../io/vtk/output_vtk.cuh"
#include "../io/debug/debug_domain.cuh"
#include "../io/diagnostics/output_tke_bin.cuh"
#include "../io/diagnostics/output_centerline_bin.cuh"
#include "../io/checkpoint/output_checkpoint.cuh"

#include "../lbm/state/lbm_state.cuh"
#include "../lbm/lbm_init_state.cuh"
#include "../lbm/lbm_mom_step.cuh"
#include "../lbm/domain/build_tags.cuh"

#include "../core/cuda_utils.cuh"
#include "../core/simulation_config.h"
#include "../core/geometry.h"

#include <chrono>
#include <iostream>
#include <algorithm>

namespace app
{
    namespace
    {
        bool should_sample_flow_averages(int step)
        {
            return step >= AVG_START_STEP;
        }

        double compute_tke_and_sample_flow_averages(const LBMState &state,
                                                    io::FlowAverages &flow_averages,
                                                    int step)
        {
            if (should_sample_flow_averages(step))
                return io::compute_ke_and_sample_flow_averages_host_2d(state, flow_averages);

            return io::compute_ke_host_2d(state);
        }
    }

    void run(const CudaConfig &cfg, const RunContext &ctx)
    {
        auto state = lbm_allocate_state();
        io::FlowAverages flow_averages = io::make_flow_averages(state.N);
        std::int64_t checkpoint_step = -1;
        int current_step = 0;
        int last_checkpoint_step = -1;

        if (ctx.restart_from_checkpoint)
        {
            const io::CheckpointConfig checkpoint_cfg = io::read_checkpoint_current(state, ctx.checkpoint_dir);
            checkpoint_step = checkpoint_cfg.step;
            current_step = static_cast<int>(std::min<std::int64_t>(checkpoint_step, N_STEPS));

            if (ctx.enable_io)
            {
                io::seed_tke_history_from_checkpoint(ctx.checkpoint_dir, ctx.out_dir,
                                                     static_cast<int>(checkpoint_step));
                if (should_sample_flow_averages(static_cast<int>(checkpoint_step)))
                    io::seed_flow_averages_from_checkpoint(ctx.checkpoint_dir, flow_averages,
                                                           static_cast<int>(checkpoint_step));
            }
        }
        else
        {
            init_state(state, cfg);
        }

        DomainTags tags = domain_tags_allocate();
        build_tags(tags);

        if (ctx.verbose)
        {
            if (ctx.show_progress)
                progress::ProgressUI::suspend_for_log();
            io::debug_domain(tags);
        }

        // ---------------- warmup ----------------
        if (!ctx.restart_from_checkpoint && ctx.enable_io)
        {
            upload_state_to_host(state);
            const double ke = compute_tke_and_sample_flow_averages(state,
                                                                   flow_averages,
                                                                   current_step);

            io::tke_bin_append(ctx.out_dir, current_step, ke);
            io::write_vti(state, cfg, current_step, ctx.out_dir);
        }

        const int t_end = N_STEPS;
        const int warmup_end = std::min(t_end, current_step + std::max(0, ctx.warmup_steps));
        while (current_step < warmup_end)
        {
            lbm_mom_step(state, cfg, tags);
            state.cur ^= 1;
            ++current_step;
        }
        CUDA_CHECK(cudaDeviceSynchronize());

        // ---------------- timers ----------------
        using clock = std::chrono::steady_clock;
        const auto wall0 = clock::now();

        const int t_begin = current_step;
        const int vti_interval = (ctx.vti_interval > 0) ? ctx.vti_interval : VTI_SAVE_INTERVAL;

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

        // ---------------- main loop ----------------
        while (current_step < t_end)
        {
            lbm_mom_step(state, cfg, tags);
            state.cur ^= 1;
            ++current_step;

            const bool save_tke = (current_step % SAVE_INTERVAL == 0);
            const bool save_vti = (current_step % vti_interval == 0);

            if (ctx.enable_io && (save_tke || save_vti))
            {
                upload_state_to_host(state);

                if (save_tke)
                {
                    const double ke = compute_tke_and_sample_flow_averages(state,
                                                                           flow_averages,
                                                                           current_step);
                    io::tke_bin_append(ctx.out_dir, current_step, ke);
                }

                if (save_vti)
                {
                    io::write_vti(state, cfg, current_step, ctx.out_dir);
                    io::write_checkpoint_current(state, current_step, ctx.out_dir);
                    if (flow_averages.samples > 0)
                        io::write_flow_averages_checkpoint(flow_averages, current_step, ctx.out_dir);
                    last_checkpoint_step = current_step;
                }
            }

            if (ctx.show_progress && (ui.should_print() || current_step == t_end))
            {
                const auto now = clock::now();

                CUDA_CHECK(cudaEventRecord(ev_prog1));
                CUDA_CHECK(cudaEventSynchronize(ev_prog1));

                float ms = 0.0f;
                CUDA_CHECK(cudaEventElapsedTime(&ms, ev_prog0, ev_prog1));
                const double gpu_elapsed_s = double(ms) * 1e-3;

                const double wall_elapsed_s = std::chrono::duration<double>(now - wall0).count();

                const int done_steps = (current_step - t_begin);
                const double updates = double(NX) * double(NY) * double(done_steps);
                const double mlups_partial = (gpu_elapsed_s > 0.0) ? (updates / gpu_elapsed_s / 1e6) : 0.0;

                ui.print(current_step, wall_elapsed_s, gpu_elapsed_s, mlups_partial);
            }
        }

        if (ctx.enable_io)
        {
            upload_state_to_host(state);
            if (last_checkpoint_step != t_end)
            {
                io::write_checkpoint_current(state, t_end, ctx.out_dir);
                if (flow_averages.samples > 0)
                    io::write_flow_averages_checkpoint(flow_averages, t_end, ctx.out_dir);
            }
            if (flow_averages.samples > 0)
                io::write_flow_averages(flow_averages, t_end, ctx.out_dir);
            io::write_centerline_profiles(state, t_end * U_LID / NX, ctx.out_dir);
        }

        const double gpu_s = gt.stop_seconds();
        const auto wall1 = clock::now();
        const double wall_s = std::chrono::duration<double>(wall1 - wall0).count();

        // ---------------- metrics ----------------
        BenchmarkResult r;
        r.gpu_seconds = gpu_s;
        r.wall_seconds = wall_s;
        r.measured_steps = (t_end - t_begin);

        const double updates = double(NX) * double(NY) * double(r.measured_steps);
        r.mlups_gpu = (r.gpu_seconds > 0.0) ? (updates / r.gpu_seconds / 1e6) : 0.0;
        r.mlups_wall = (r.wall_seconds > 0.0) ? (updates / r.wall_seconds / 1e6) : 0.0;

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
