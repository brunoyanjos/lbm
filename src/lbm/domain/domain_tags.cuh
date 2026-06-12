#pragma once
#include <cstdint>
#include <cstddef>

#include "core/mask_type.cuh"

enum class NodeId : uint8_t
{
    SOLID = 0,
    NORTH = 3,
    NORTH_EAST = 1,
    NORTH_WEST = 2,
    EAST = 5,
    WEST = 10,
    SOUTH = 12,
    SOUTH_EAST = 4,
    SOUTH_WEST = 8,
    FLUID = 15,
};

struct DomainTags
{
    // device
    mask_t *d_valid = nullptr;
    uint8_t *d_node = nullptr;

    // host
    mask_t *h_valid = nullptr;
    uint8_t *h_node = nullptr;

    size_t N = 0;
    size_t bytes_valid = 0;
    size_t bytes_node = 0;
};

DomainTags domain_tags_allocate();
void domain_tags_free(DomainTags &T);
