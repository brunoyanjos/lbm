#pragma once

#include "lbm_features.cuh"
#include "../lbm/stencil_active.cuh"

namespace lbm
{
    static constexpr int reg_order = LBM_REG_ORDER;

    static constexpr bool stencil_is_d2q9 = (LBM_ACTIVE_STENCIL_D2Q9 != 0);
    static constexpr bool stencil_is_d2v17 = (LBM_ACTIVE_STENCIL_D2V17 != 0);

    static constexpr bool has_reg3_cross = (LBM_HAS_REG3_CROSS != 0);
    static constexpr bool has_reg3_axial = (LBM_HAS_REG3_AXIAL != 0);
}