#pragma once

#include "../../config/reg_order.cuh"
#include "../../stencils/stencil_traits.cuh"

namespace lbm_reg
{

    template <int RegOrder, class StencilTag>
    struct Traits
    {
        static constexpr bool has_reg3_cross = false;
        static constexpr bool has_reg3_axial = false;
    };

    template <int RegOrder>
    struct Traits<RegOrder, D2Q9Tag>
    {
        static constexpr bool has_reg3_cross =
            (RegOrder >= 3);

        static constexpr bool has_reg3_axial =
            false;
    };

    template <int RegOrder>
    struct Traits<RegOrder, D2V17Tag>
    {
        static constexpr bool has_reg3_cross =
            (RegOrder >= 3);

        static constexpr bool has_reg3_axial =
            (RegOrder >= 3);
    };

    using Active = Traits<
        lbm_config::reg_order,
        ActiveStencilTag>;

}