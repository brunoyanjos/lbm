#pragma once
#include <cstdint>
#include <cstddef>

enum class WallId : uint8_t
{
    NONE = 0,
    LEFT = 1,
    RIGHT = 2,
    BOTTOM = 3,
    TOP = 4,
    CORNER = 5
};

struct DomainTags
{
    // device
    uint32_t *d_valid = nullptr; // bit i = 1 => direção i é válida (vizinho dentro do domínio)
    uint8_t *d_wall = nullptr;   // WallId por nó (NONE para bulk)

    // host optional (debug)
    uint32_t *h_valid = nullptr;
    uint8_t *h_wall = nullptr;

    size_t N = 0;
    size_t bytes_valid = 0;
    size_t bytes_wall = 0;
};

DomainTags domain_tags_allocate(bool host_buffers);
void domain_tags_free(DomainTags &T);
