#pragma once

#if defined(LBM_GEOM_ANNUL)
#include "annul/build_tags.cuh"
#include "annul/bc_velocity.cuh"
#include "annul/outputs.cuh"
#include "annul/properties.cuh"

namespace Geometry = ANNUL;

#elif defined(LBM_GEOM_CHANNEL)
#include "poiseuille/build_tags.cuh"
#include "poiseuille/bc_velocity.cuh"
#include "poiseuille/outputs.cuh"
#include "poiseuille/properties.cuh"

namespace Geometry = POISEUILLE;

#elif defined(LBM_GEOM_COUETTE)
#include "couette/build_tags.cuh"
#include "couette/bc_velocity.cuh"
#include "couette/outputs.cuh"
#include "couette/properties.cuh"

namespace Geometry = COUETTE;

#elif defined(LBM_GEOM_JET)
#include "jet/build_tags.cuh"
#include "jet/bc_velocity.cuh"
#include "jet/outputs.cuh"
#include "jet/properties.cuh"

namespace Geometry = JET;

#elif defined(LBM_GEOM_SQUARE_CAVITY)
#include "square_cavity/build_tags.cuh"
#include "square_cavity/bc_velocity.cuh"
#include "square_cavity/outputs.cuh"
#include "square_cavity/properties.cuh"

namespace Geometry = SQUARE_CAVITY;

#else
#include "square_cavity/build_tags.cuh"
#include "square_cavity/bc_velocity.cuh"
#include "square_cavity/outputs.cuh"
#include "square_cavity/properties.cuh"

namespace Geometry = SQUARE_CAVITY;

#endif
