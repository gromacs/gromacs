#include "gromacs/utility/basedefinitions.h"

#define R123_CUDA_DEVICE /* Random123 isn't used from Cuda */
#define R123_STATIC_INLINE static gmx_inline
/* force_inline isn't used in gromacs - if it matters for a compiler it probably
    not only matters here and should be defined in basedefinitions */
#define R123_FORCE_INLINE(decl) decl 
#include <assert.h>
#define R123_ASSERT assert
#define R123_STATIC_ASSERT(expr, msg) /* Not used by the code we use */
