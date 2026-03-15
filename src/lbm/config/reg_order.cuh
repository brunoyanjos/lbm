#pragma once

#ifndef LBM_REG_ORDER
#define LBM_REG_ORDER 2
#endif

static_assert(LBM_REG_ORDER == 2 || LBM_REG_ORDER == 3,
              "LBM_REG_ORDER must be 2 or 3");

namespace lbm_config
{
    inline constexpr int reg_order = LBM_REG_ORDER;
}