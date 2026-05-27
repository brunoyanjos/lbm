#pragma once

#include <cstddef>

template <auto... Ids>
struct IdList
{
    static constexpr std::size_t size = sizeof...(Ids);
};