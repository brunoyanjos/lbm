#pragma once

#include "lbm/boundary/common/id_list.cuh"
#include "lbm/moment/moment_id.cuh"

using DirichletUnknownMoments = IdList<
    MomentId::mxy>;
