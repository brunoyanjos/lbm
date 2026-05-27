#pragma once

#include "lbm/boundary/common/id_list.cuh"

using NonLinearSystemList = IdList<
    NonlinearMomentId::uxux,
    NonlinearMomentId::uxuy,
    NonlinearMomentId::uyuy,
    >;
