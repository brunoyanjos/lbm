#pragma once
#include "../lbm/domain/domain_tags.cuh"

namespace io
{
    void debug_domain(const DomainTags &T,
                      bool print_domain = true,
                      bool print_masks = true,
                      int max_mask_points = 2000);
}