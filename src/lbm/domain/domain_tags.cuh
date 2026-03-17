#pragma once
#include <cstdint>
#include <cstddef>

enum class NodeId : uint8_t
{
    FLUID = 0,
    SOLID = 1,
    DIRICHLET = 2,
    INLET = 3,
    OUTLET = 4,
};

struct DomainTags
{
    // device
    uint32_t *d_valid = nullptr;
    uint8_t *d_node = nullptr;

    // host
    uint32_t *h_valid = nullptr;
    uint8_t *h_node = nullptr;

    size_t N = 0;
    size_t bytes_valid = 0;
    size_t bytes_node = 0;
};

DomainTags domain_tags_allocate();
void domain_tags_free(DomainTags &T);
