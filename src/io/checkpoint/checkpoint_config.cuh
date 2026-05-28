#pragma once

#include "../../core/geometry.h"
#include "../../core/physics.h"
#include "../../core/simulation_config.h"
#include "../../core/types.cuh"
#include "../../lbm/stencil_active.cuh"

#include <cstdint>
#include <cstring>
#include <type_traits>

namespace io
{
    enum class CheckpointStencil : std::uint32_t
    {
        D2Q9 = 9,
        D2V17 = 17,
        D2V37 = 37,
    };

    struct CheckpointConfig
    {
        char magic[8];
        std::uint32_t version;
        std::uint32_t header_bytes;
        std::uint32_t endian_marker;

        std::int32_t nx;
        std::int32_t ny;
        std::int32_t q;
        std::int32_t real_bytes;
        std::int32_t real_is_double;
        std::int32_t cur;
        std::int32_t field_count;
        std::int32_t reg_order;
        std::int32_t recurrence;

        std::int64_t node_count;
        std::int64_t field_bytes;
        std::int64_t payload_bytes;
        std::int64_t step;
        std::int64_t n_steps;
        std::int64_t save_interval;
        std::int64_t vti_save_interval;

        double re;
        double u_lid;
        double tau;
        double omega;

        CheckpointStencil stencil;
        char stencil_name[16];
    };

    static_assert(std::is_trivially_copyable_v<CheckpointConfig>,
                  "CheckpointConfig must stay binary-serializable");

    inline CheckpointStencil active_checkpoint_stencil()
    {
#if defined(LBM_STENCIL_D2V37)
        return CheckpointStencil::D2V37;
#elif defined(LBM_STENCIL_D2V17)
        return CheckpointStencil::D2V17;
#else
        return CheckpointStencil::D2Q9;
#endif
    }

    inline const char *active_checkpoint_stencil_name()
    {
#if defined(LBM_STENCIL_D2V37)
        return "D2V37";
#elif defined(LBM_STENCIL_D2V17)
        return "D2V17";
#else
        return "D2Q9";
#endif
    }

    inline CheckpointConfig make_checkpoint_config(std::int64_t step, int cur)
    {
        CheckpointConfig cfg{};
        std::memcpy(cfg.magic, "LBMCHK1", 8);
        cfg.version = 2;
        cfg.header_bytes = sizeof(CheckpointConfig);
        cfg.endian_marker = 0x01020304u;

        cfg.nx = NX;
        cfg.ny = NY;
        cfg.q = Stencil::Q;
        cfg.real_bytes = sizeof(real_t);
        cfg.real_is_double = (sizeof(real_t) == sizeof(double)) ? 1 : 0;
        cfg.cur = cur;
        cfg.reg_order = REG_ORDER;
        cfg.recurrence = USE_RECURRENCE ? 1 : 0;
        cfg.field_count = 6;
        if constexpr (REG_ORDER == 3 && !USE_RECURRENCE)
            cfg.field_count += Stencil::high_order ? 4 : 2;

        cfg.node_count = static_cast<std::int64_t>(NX) * static_cast<std::int64_t>(NY);
        cfg.field_bytes = cfg.node_count * static_cast<std::int64_t>(sizeof(real_t));
        cfg.payload_bytes = cfg.field_count * cfg.field_bytes;
        cfg.step = step;
        cfg.n_steps = N_STEPS;
        cfg.save_interval = SAVE_INTERVAL;
        cfg.vti_save_interval = VTI_SAVE_INTERVAL;

        cfg.re = static_cast<double>(RE);
        cfg.u_lid = static_cast<double>(U_LID);
        cfg.tau = static_cast<double>(TAU);
        cfg.omega = static_cast<double>(OMEGA);

        cfg.stencil = active_checkpoint_stencil();
        std::strncpy(cfg.stencil_name, active_checkpoint_stencil_name(),
                     sizeof(cfg.stencil_name) - 1);

        return cfg;
    }
} // namespace io
