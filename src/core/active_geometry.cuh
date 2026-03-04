#pragma once

#if defined(__INTELLISENSE__) &&                              \
    !defined(LBM_GEOM_ANNUL) && !defined(LBM_GEOM_CHANNEL) && \
    !defined(LBM_GEOM_COUETTE) && !defined(LBM_GEOM_JET) &&   \
    !defined(LBM_GEOM_SQUARE_CAVITY)
#define LBM_GEOM_SQUARE_CAVITY 1
#endif

#if (defined(LBM_GEOM_ANNUL) + defined(LBM_GEOM_CHANNEL) + defined(LBM_GEOM_COUETTE) + \
     defined(LBM_GEOM_JET) + defined(LBM_GEOM_SQUARE_CAVITY)) != 1
#error "Select exactly one geometry: -DLBM_GEOM_ANNUL / CHANNEL / COUETTE / JET / SQUARE_CAVITY"
#endif

#if defined(LBM_GEOM_ANNUL)
#include "geometries/annul/geometry.h"
#include "geometries/annul/physics.h"

#elif defined(LBM_GEOM_CHANNEL)
#include "geometries/poiseuille/geometry.h"
#include "geometries/poiseuille/physics.h"

#elif defined(LBM_GEOM_COUETTE)
#include "geometries/couette/geometry.h"
#include "geometries/couette/physics.h"

#elif defined(LBM_GEOM_JET)
#include "geometries/jet/geometry.h"
#include "geometries/jet/physics.h"

#elif defined(LBM_GEOM_SQUARE_CAVITY)
#include "geometries/square_cavity/geometry.h"
#include "geometries/square_cavity/physics.h"

#endif