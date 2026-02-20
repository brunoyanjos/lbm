#pragma once
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace progress
{
    struct ProgressUI
    {
        int t0 = 0;
        int t_end = 0;

        bool enabled = false;
        double hz = 2.0;

        std::chrono::steady_clock::time_point wall0{};
        std::chrono::steady_clock::time_point last_print{};

        // eventos CUDA opcionais (você já mede GPU parcial no simulation)
        // aqui a UI só imprime

        void start(int t0_, int t_end_, bool enabled_, double hz_)
        {
            t0 = t0_;
            t_end = t_end_;
            enabled = enabled_;
            hz = (hz_ > 0.0 ? hz_ : 2.0);

            wall0 = std::chrono::steady_clock::now();
            last_print = wall0;
        }

        static inline std::string fmt_hms(double s)
        {
            int si = (s < 0) ? 0 : int(s + 0.5);
            int h = si / 3600;
            si %= 3600;
            int m = si / 60;
            si %= 60;
            int sec = si;

            std::ostringstream os;
            os << h << ":" << std::setw(2) << std::setfill('0') << m
               << ":" << std::setw(2) << std::setfill('0') << sec;
            return os.str();
        }

        static inline void clear_line()
        {
            // limpa a linha atual no terminal (ANSI)
            // \r -> início da linha
            // \033[2K -> clear entire line
            std::cerr << "\r\033[2K" << std::flush;
        }

        // chame antes de imprimir qualquer log em stdout (ou stderr)
        static inline void suspend_for_log()
        {
            // encerra a linha do progresso pra não misturar
            std::cerr << "\r\033[2K" << std::flush;
        }

        bool should_print() const
        {
            if (!enabled)
                return false;
            const auto now = std::chrono::steady_clock::now();
            const double dt = std::chrono::duration<double>(now - last_print).count();
            const double min_dt = (hz > 0.0) ? (1.0 / hz) : 0.5;
            return dt >= min_dt;
        }

        void print(int t,
                   double wall_elapsed_s,
                   double gpu_elapsed_s,
                   double mlups_partial)
        {
            if (!enabled)
                return;

            const int total = (t_end - t0);
            const int done = std::clamp(t - t0, 0, total);
            const double frac = (total > 0) ? (double(done) / double(total)) : 1.0;
            const double eta = (frac > 1e-12) ? (wall_elapsed_s * (1.0 / frac - 1.0)) : 0.0;

            const int bar_w = 40;
            const int pos = int(bar_w * frac);

            std::ostringstream os;
            os << "\r[";
            for (int i = 0; i < bar_w; ++i)
                os << (i < pos ? "=" : " ");
            os << "] ";

            os << std::setw(3) << int(frac * 100.0) << "% "
               << "step " << t << "/" << t_end
               << " | wall " << fmt_hms(wall_elapsed_s)
               << " | ETA " << fmt_hms(eta)
               << " | gpu " << std::fixed << std::setprecision(3) << gpu_elapsed_s << "s"
               << " | MLUPS " << std::fixed << std::setprecision(2) << mlups_partial
               << std::flush;

            std::cerr << os.str();

            last_print = std::chrono::steady_clock::now();
        }

        // limpa e não deixa “resto” da barra na tela
        void finish(bool keep_final_line = false,
                    const char *final_msg = nullptr)
        {
            if (!enabled)
                return;

            if (keep_final_line)
            {
                // imprime uma linha final (sem \r) e pula linha
                clear_line();
                if (final_msg)
                    std::cerr << final_msg;
                std::cerr << "\n";
            }
            else
            {
                // “some” com o progresso
                clear_line();
            }
        }
    };
} // namespace progress
