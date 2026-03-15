#pragma once

struct D2Q9Tag
{
};

struct D2V17Tag
{
};

#if defined(LBM_STENCIL_D2Q9)
using ActiveStencilTag = D2Q9Tag;
#elif defined(LBM_STENCIL_D2V17)
using ActiveStencilTag = D2V17Tag;
#else
using ActiveStencilTag = D2Q9Tag;
#endif