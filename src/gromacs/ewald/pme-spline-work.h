#ifndef GMX_EWALD_PME_SPLINE_WORK
#define GMX_EWALD_PME_SPLINE_WORK

#ifdef __cplusplus
extern "C" {
#endif

#include "pme-simd.h"

struct pme_spline_work
{
#ifdef PME_SIMD4_SPREAD_GATHER
    /* Masks for 4-wide SIMD aligned spreading and gathering */
    gmx_simd4_bool_t mask_S0[6], mask_S1[6];
#else
    int              dummy; /* C89 requires that struct has at least one member */
#endif
};

struct pme_spline_work *make_pme_spline_work(int gmx_unused order);

#ifdef __cplusplus
}
#endif

#endif

