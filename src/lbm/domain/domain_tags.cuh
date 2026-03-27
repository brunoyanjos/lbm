#pragma once
#include <cstdint>
#include <cstddef>

#include "core/mask_type.cuh"

enum class NodeId : uint8_t
{
    SOLID = 0,
    ONE = 1,
    TWO = 2,
    THREE = 3,
    FOUR = 4,
    FIVE = 5,
    SEVEN = 7,
    EIGHT = 8,
    TEN = 10,
    ELEVEN = 11,
    TWELVE = 12,
    THIRTEEN = 13,
    FOURTEEN = 14,
    FLUID = 15
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
