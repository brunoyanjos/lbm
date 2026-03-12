#pragma once

#include "id_list.cuh"

template <auto Query, class List>
struct IndexOfId;

template <auto Query>
struct IndexOfId<Query, IdList<>>
{
    static constexpr int value = -1;
};

template <auto Query, auto First, auto... Rest>
struct IndexOfId<Query, IdList<First, Rest...>>
{
private:
    static constexpr int tail_value = IndexOfId<Query, IdList<Rest...>>::value;

public:
    static constexpr int value =
        (Query == First) ? 0 : (tail_value < 0 ? -1 : 1 + tail_value);
};

template <auto Query, class List>
inline constexpr int index_of_id_v = IndexOfId<Query, List>::value;

template <auto Query, class List>
inline constexpr bool has_id_v = (index_of_id_v<Query, List> >= 0);