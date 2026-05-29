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
#include <vector>

namespace app
{
    namespace
    {
        bool should_sample_flow_averages(int step)
        {
            return step >= AVG_START_STEP;
        }

        double compute_tke_and_sample_flow_averages(const std::vector<LBMState> &state,
                                                    io::FlowAverages &flow_averages,
                                                    int step)
        {
            if (should_sample_flow_averages(step))
                return io::compute_ke_and_sample_flow_averages_host_2d(state, flow_averages);

            return io::compute_ke_host_2d(state);
        }

        void copy_device_rows(real_t *dst, int dst_device, int dst_y,
                              const real_t *src, int src_device, int src_y,
                              int rows, int nx)
        {
            const size_t dst_offset = static_cast<size_t>(dst_y) * static_cast<size_t>(nx);
            const size_t src_offset = static_cast<size_t>(src_y) * static_cast<size_t>(nx);
            const size_t bytes = static_cast<size_t>(rows) * static_cast<size_t>(nx) * sizeof(real_t);

            if (dst_device == src_device)
            {
                CUDA_CHECK(cudaSetDevice(dst_device));
                CUDA_CHECK(cudaMemcpy(dst + dst_offset,
                                      src + src_offset,
                                      bytes,
                                      cudaMemcpyDeviceToDevice));
            }
            else
            {
                CUDA_CHECK(cudaMemcpyPeer(dst + dst_offset,
                                          dst_device,
                                          src + src_offset,
                                          src_device,
                                          bytes));
            }
        }

        template <int RegOrder, bool Rec, bool HighOrder>
        void exchange_extra_halo_fields(LBMStateFor<RegOrder, Rec, HighOrder> &,
                                        int,
                                        int,
                                        const LBMStateFor<RegOrder, Rec, HighOrder> &,
                                        int,
                                        const LBMStateFor<RegOrder, Rec, HighOrder> &,
                                        int)
        {
        }

        template <bool HighOrder>
        void exchange_extra_halo_fields(LBMStateFor<3, false, HighOrder> &state,
                                        int device,
                                        int buffer,
                                        const LBMStateFor<3, false, HighOrder> &prev,
                                        int prev_device,
                                        const LBMStateFor<3, false, HighOrder> &next,
                                        int next_device)
        {
            const int halo = state.domain.halo;
            const int nx = state.domain.nx;
            const int lower_dst = 0;
            const int upper_dst = state.domain.halo + state.domain.local_ny;
            const int prev_src = prev.domain.halo + prev.domain.local_ny - halo;
            const int next_src = next.domain.halo;

            copy_device_rows(state.d_mxxy[buffer], device, lower_dst, prev.d_mxxy[buffer], prev_device, prev_src, halo, nx);
            copy_device_rows(state.d_mxxy[buffer], device, upper_dst, next.d_mxxy[buffer], next_device, next_src, halo, nx);
            copy_device_rows(state.d_mxyy[buffer], device, lower_dst, prev.d_mxyy[buffer], prev_device, prev_src, halo, nx);
            copy_device_rows(state.d_mxyy[buffer], device, upper_dst, next.d_mxyy[buffer], next_device, next_src, halo, nx);

            if constexpr (HighOrder)
            {
                copy_device_rows(state.d_mxxx[buffer], device, lower_dst, prev.d_mxxx[buffer], prev_device, prev_src, halo, nx);
                copy_device_rows(state.d_mxxx[buffer], device, upper_dst, next.d_mxxx[buffer], next_device, next_src, halo, nx);
                copy_device_rows(state.d_myyy[buffer], device, lower_dst, prev.d_myyy[buffer], prev_device, prev_src, halo, nx);
                copy_device_rows(state.d_myyy[buffer], device, upper_dst, next.d_myyy[buffer], next_device, next_src, halo, nx);
            }
        }

        void exchange_halos(std::vector<LBMState> &states, const RunContext &ctx)
        {
            const size_t n = states.size();
            if (n <= 1)
                return;

            for (size_t i = 0; i < n; ++i)
            {
                LBMState &state = states[i];
                const LBMState &prev = states[(i + n - 1) % n];
                const LBMState &next = states[(i + 1) % n];

                const int device = ctx.partitions[i].device_id;
                const int prev_device = ctx.partitions[(i + n - 1) % n].device_id;
                const int next_device = ctx.partitions[(i + 1) % n].device_id;
                const int buffer = state.cur;
                const int halo = state.domain.halo;
                const int nx = state.domain.nx;
                const int lower_dst = 0;
                const int upper_dst = state.domain.halo + state.domain.local_ny;
                const int prev_src = prev.domain.halo + prev.domain.local_ny - halo;
                const int next_src = next.domain.halo;

                copy_device_rows(state.d_rho[buffer], device, lower_dst, prev.d_rho[buffer], prev_device, prev_src, halo, nx);
                copy_device_rows(state.d_rho[buffer], device, upper_dst, next.d_rho[buffer], next_device, next_src, halo, nx);
                copy_device_rows(state.d_ux[buffer], device, lower_dst, prev.d_ux[buffer], prev_device, prev_src, halo, nx);
                copy_device_rows(state.d_ux[buffer], device, upper_dst, next.d_ux[buffer], next_device, next_src, halo, nx);
                copy_device_rows(state.d_uy[buffer], device, lower_dst, prev.d_uy[buffer], prev_device, prev_src, halo, nx);
                copy_device_rows(state.d_uy[buffer], device, upper_dst, next.d_uy[buffer], next_device, next_src, halo, nx);
                copy_device_rows(state.d_mxx[buffer], device, lower_dst, prev.d_mxx[buffer], prev_device, prev_src, halo, nx);
                copy_device_rows(state.d_mxx[buffer], device, upper_dst, next.d_mxx[buffer], next_device, next_src, halo, nx);
                copy_device_rows(state.d_mxy[buffer], device, lower_dst, prev.d_mxy[buffer], prev_device, prev_src, halo, nx);
                copy_device_rows(state.d_mxy[buffer], device, upper_dst, next.d_mxy[buffer], next_device, next_src, halo, nx);
                copy_device_rows(state.d_myy[buffer], device, lower_dst, prev.d_myy[buffer], prev_device, prev_src, halo, nx);
                copy_device_rows(state.d_myy[buffer], device, upper_dst, next.d_myy[buffer], next_device, next_src, halo, nx);

                exchange_extra_halo_fields(state, device, buffer, prev, prev_device, next, next_device);
            }
        }
    }

    void run(const RunContext &ctx)
    {
        auto states = allocate_partition_states(ctx);

        io::FlowAverages flow_averages = io::make_flow_averages(static_cast<size_t>(NX) * static_cast<size_t>(NY));
        std::int64_t checkpoint_step = -1;
        int current_step = 0;
        int last_checkpoint_step = -1;

        // if (ctx.restart_from_checkpoint)
        // {
        //     const io::CheckpointConfig checkpoint_cfg = io::read_checkpoint_current(state, ctx.checkpoint_dir);
        //     checkpoint_step = checkpoint_cfg.step;
        //     current_step = static_cast<int>(std::min<std::int64_t>(checkpoint_step, N_STEPS));

        //     if (ctx.enable_io)
        //     {
        //         io::seed_tke_history_from_checkpoint(ctx.checkpoint_dir, ctx.out_dir,
        //                                              static_cast<int>(checkpoint_step));
        //         if (should_sample_flow_averages(static_cast<int>(checkpoint_step)))
        //             io::seed_flow_averages_from_checkpoint(ctx.checkpoint_dir, flow_averages,
        //                                                    static_cast<int>(checkpoint_step));
        //     }
        // }
        // else
        // {
        init_state(states, ctx);
        exchange_halos(states, ctx);
        // }

        auto tags = allocate_partition_tags(ctx);

        // if (ctx.verbose)
        // {
        //     if (ctx.show_progress)
        //         progress::ProgressUI::suspend_for_log();
        //     io::debug_domain(tags);
        // }

        // // ---------------- warmup ----------------
        // if (!ctx.restart_from_checkpoint && ctx.enable_io)
        // {
        upload_state_to_host(states, ctx);
        const double ke = compute_tke_and_sample_flow_averages(states,
                                                               flow_averages,
                                                               current_step);

        io::tke_bin_append(ctx.out_dir, current_step, ke);
        // }

        const int t_end = N_STEPS;
        const int warmup_end = std::min(t_end, current_step + std::max(0, ctx.warmup_steps));
        // while (current_step < warmup_end)
        // {
        //     lbm_mom_step(state, cfg, tags);
        //     state.cur ^= 1;
        //     ++current_step;
        // }
        // CUDA_CHECK(cudaDeviceSynchronize());

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
            exchange_halos(states, ctx);
            lbm_mom_step(states, tags, ctx);
            ++current_step;

            const bool save_tke = (current_step % SAVE_INTERVAL == 0);
            const bool save_vti = (current_step % vti_interval == 0);

            if (ctx.enable_io && (save_tke || save_vti))
            {
                upload_state_to_host(states, ctx);

                if (save_tke)
                {
                    const double ke = compute_tke_and_sample_flow_averages(states,
                                                                           flow_averages,
                                                                           current_step);
                    io::tke_bin_append(ctx.out_dir, current_step, ke);
                }

                if (save_vti)
                {
                    io::write_vti(states, current_step, ctx.out_dir);
                    // io::write_checkpoint_current(state, current_step, ctx.out_dir);
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
            upload_state_to_host(states, ctx);
            if (last_checkpoint_step != t_end)
            {
                // io::write_checkpoint_current(state, t_end, ctx.out_dir);
                if (flow_averages.samples > 0)
                    io::write_flow_averages_checkpoint(flow_averages, t_end, ctx.out_dir);
            }
            if (flow_averages.samples > 0)
                io::write_flow_averages(flow_averages, t_end, ctx.out_dir);
            // io::write_centerline_profiles(state, t_end * U_LID / NX, ctx.out_dir);
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

        io::write_performance(ctx.out_dir, r);

        ui.finish(false);

        if (ctx.show_progress)
        {
            CUDA_CHECK(cudaEventDestroy(ev_prog0));
            CUDA_CHECK(cudaEventDestroy(ev_prog1));
        }

        std::cout << "Simulation finished after " << N_STEPS << " timesteps.\n";
        std::cout << "GPU(s)=" << r.gpu_seconds << " WALL(s)=" << r.wall_seconds
                  << " MLUPS_GPU=" << r.mlups_gpu << " MLUPS_WALL=" << r.mlups_wall << "\n";

        // domain_tags_free(tags);
        free_partition_tags(tags, ctx);
        free_partition_states(states, ctx);
    }
} // namespace app
