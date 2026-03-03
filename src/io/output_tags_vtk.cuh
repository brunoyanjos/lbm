#pragma once
#include <string>
#include "../lbm/domain/active_geometry.cuh"

namespace io
{
    void write_tags_vtk(const DomainTags &T, const std::string &out_dir);
}
