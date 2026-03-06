#pragma once

#if defined(__INTELLISENSE__) &&                              \
    !defined(LBM_GEOM_ANNUL) && !defined(LBM_GEOM_CHANNEL) && \
    !defined(LBM_GEOM_COUETTE) && !defined(LBM_GEOM_JET) &&   \
    !defined(LBM_GEOM_SQUARE_CAVITY)
#define LBM_GEOM_ANNUL 1
#endif

#if (defined(LBM_GEOM_ANNUL) + defined(LBM_GEOM_CHANNEL) + defined(LBM_GEOM_COUETTE) + \
     defined(LBM_GEOM_JET) + defined(LBM_GEOM_SQUARE_CAVITY)) != 1
#error "Select exactly one geometry: -DLBM_GEOM_ANNUL / CHANNEL / COUETTE / JET / SQUARE_CAVITY"
#endif

#if defined(LBM_GEOM_ANNUL)
#include "geometries/annul/build_tags.cuh"
#include "geometries/annul/bc_velocity.cuh"
namespace Geometry = ANNUL;

#elif defined(LBM_GEOM_CHANNEL)
#include "geometries/poiseuille/build_tags.cuh"
#include "geometries/poiseuille/bc_velocity.cuh"
namespace Geometry = POISEUILLE;

#elif defined(LBM_GEOM_COUETTE)
#include "geometries/couette/build_tags.cuh"
#include "geometries/couette/bc_velocity.cuh"
namespace Geometry = COUETTE;

#elif defined(LBM_GEOM_JET)
#include "geometries/jet/build_tags.cuh"
#include "geometries/jet/bc_velocity.cuh"
namespace Geometry = JET;

#elif defined(LBM_GEOM_SQUARE_CAVITY)
#include "geometries/square_cavity/build_tags.cuh"
#include "geometries/square_cavity/bc_velocity.cuh"
namespace Geometry = SQUARE_CAVITY;

#endif