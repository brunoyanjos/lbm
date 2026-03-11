#pragma once

#include <cstdint>
#include <type_traits>

#if defined(REAL_T_IS_DOUBLE)
using real_t = double;
#else
using real_t = float;
#endif

template <class E>
__host__ __device__ __forceinline__ constexpr std::enable_if_t<std::is_enum_v<E>, uint8_t> to_u8(E e)
{
    return static_cast<uint8_t>(e);
}

template <class T>
__host__ __device__ __forceinline__ constexpr real_t r_cast(T v) { return static_cast<real_t>(v); }

namespace r
{
    static constexpr real_t zero = real_t{0};
    static constexpr real_t one = real_t{1};
    static constexpr real_t two = real_t{2};
    static constexpr real_t three = real_t{3};

    static constexpr real_t quarter = real_t{0.25};
    static constexpr real_t third = real_t{1.0 / 3.0};
    static constexpr real_t half = real_t{0.5};
    static constexpr real_t sixth = real_t{1.0 / 6.0};
}