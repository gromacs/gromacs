#include "gromacs/utility/basedefinitions.h"
#define R123_CUDA_DEVICE
#define R123_STATIC_INLINE static gmx_inline
#define R123_FORCE_INLINE(decl) decl
#include <assert.h>
#define R123_ASSERT assert
#define R123_STATIC_ASSERT(expr, msg)
