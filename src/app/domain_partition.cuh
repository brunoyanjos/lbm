#pragma once

#include "../core/geometry.h"
#include "../lbm/stencil_active.cuh"

#include <vector>

namespace app
{

    struct DomainPartition
    {
        int device_id = 0;
        int y_begin = 0;
        int y_end = NY;
        int local_ny = NY;
        int halo = Stencil::radius;
    };

    inline std::vector<DomainPartition> make_domain_partitions(const std::vector<int> &device_ids)
    {
        std::vector<DomainPartition> partitions;
        partitions.reserve(device_ids.size());

        const int count = static_cast<int>(device_ids.size());
        const int base_ny = NY / count;
        const int extra = NY % count;

        int y = 0;
        for (int i = 0; i < count; ++i)
        {
            const int local_ny = base_ny + (i < extra ? 1 : 0);

            DomainPartition p{};
            p.device_id = device_ids[static_cast<size_t>(i)];
            p.y_begin = y;
            p.y_end = y + local_ny;
            p.local_ny = local_ny;
            p.halo = Stencil::radius;
            partitions.push_back(p);

            y = p.y_end;
        }

        return partitions;
    }

    inline std::vector<DomainPartition> make_single_device_partition(int device_id)
    {
        return make_domain_partitions(std::vector<int>{device_id});
    }

} // namespace app
