#pragma once

template <auto... Ids>
struct IdList
{
    static constexpr int size = sizeof...(Ids);
};