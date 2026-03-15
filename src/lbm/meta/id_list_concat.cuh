#pragma once

#include "id_list.cuh"

template <class A, class B>
struct IdListConcat;

template <auto... AIds, auto... BIds>
struct IdListConcat<IdList<AIds...>, IdList<BIds...>>
{
    using type = IdList<AIds..., BIds...>;
};

template <class A, class B>
using IdListConcatT = typename IdListConcat<A, B>::type;