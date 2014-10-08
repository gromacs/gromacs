#ifndef GMX_EWALD_PME_GATHER
#define GMX_EWALD_PME_GATHER

#ifdef __cplusplus
extern "C" {
#endif

#include "pme-internal.h"
#include "gromacs/utility/real.h"

void
gather_f_bsplines(struct gmx_pme *pme, real *grid,
                  gmx_bool bClearF, pme_atomcomm_t *atc,
                  splinedata_t *spline,
                  real scale);

real
gather_energy_bsplines(struct gmx_pme *pme, real *grid,
                       pme_atomcomm_t *atc);

#ifdef __cplusplus
}
#endif

#endif

