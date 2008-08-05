#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef _perf_est_h
#define _perf_est_h

#include "typedefs.h"

extern int n_bonded_dx(gmx_mtop_t *mtop,bool bExcl);
/* Returns the number of pbc_rvec_sub calls required for bonded interactions.
 * This number is also roughly proportional to the computational cost.
 */

extern float pme_load_estimate(gmx_mtop_t *mtop,t_inputrec *ir,matrix box);
/* Returns an estimate for the relative load of the PME mesh calculation
 * in the total force calculation.
 * This estimate is reasonable for recent Intel and AMD x86_64 CPUs.
 */

#endif	/* _perf_est_h */
