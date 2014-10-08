#ifndef GMX_EWALD_PME_SIMD_H
#define GMX_EWALD_PME_SIMD_H

/* Include the SIMD macro file and then check for support */
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#ifdef GMX_SIMD_HAVE_REAL
/* Turn on arbitrary width SIMD intrinsics for PME solve */
#    define PME_SIMD_SOLVE
#endif

/* Check if we have 4-wide SIMD macro support */
#if (defined GMX_SIMD4_HAVE_REAL)
/* Do PME spread and gather with 4-wide SIMD.
 * NOTE: SIMD is only used with PME order 4 and 5 (which are the most common).
 */
#    define PME_SIMD4_SPREAD_GATHER

#    if (defined GMX_SIMD_HAVE_LOADU) && (defined GMX_SIMD_HAVE_STOREU)
/* With PME-order=4 on x86, unaligned load+store is slightly faster
 * than doubling all SIMD operations when using aligned load+store.
 */
#        define PME_SIMD4_UNALIGNED
#    endif
#endif

#ifdef PME_SIMD4_SPREAD_GATHER
#    define SIMD4_ALIGNMENT  (GMX_SIMD4_WIDTH*sizeof(real))
#else
/* We can use any alignment, apart from 0, so we use 4 reals */
#    define SIMD4_ALIGNMENT  (4*sizeof(real))
#endif

#endif
