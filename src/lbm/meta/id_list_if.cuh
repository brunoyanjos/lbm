#pragma once

template <bool Cond, class TrueList, class FalseList>
struct IdListIf
{
    using type = TrueList;
};

template <class TrueList, class FalseList>
struct IdListIf<false, TrueList, FalseList>
{
    using type = FalseList;
};

template <bool Cond, class TrueList, class FalseList>
using IdListIfT = typename IdListIf<Cond, TrueList, FalseList>::type;