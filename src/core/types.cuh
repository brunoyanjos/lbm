#pragma once

#if defined(REAL_T_IS_DOUBLE)
using real_t = double;
#else
using real_t = float;
#endif