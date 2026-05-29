#pragma once

#include "core/geometry.h"

#include <cuda_runtime.h>
#include <cstddef>

struct LocalDomain
{
    int nx = NX;
    int ny_global = NY;
    int y_begin = 0;
    int y_end = NY;
    int local_ny = NY;
    int halo = 0;
    int storage_ny = NY;
    size_t N = static_cast<size_t>(NX) * static_cast<size_t>(NY);
};

inline LocalDomain make_local_domain(int y_begin, int y_end, int halo)
{
    LocalDomain d{};
    d.nx = NX;
    d.ny_global = NY;
    d.y_begin = y_begin;
    d.y_end = y_end;
    d.local_ny = y_end - y_begin;
    d.halo = halo;
    d.storage_ny = d.local_ny + 2 * d.halo;
    d.N = static_cast<size_t>(d.nx) * static_cast<size_t>(d.storage_ny);
    return d;
}
