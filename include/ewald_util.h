

#ifndef _ewald_util_h
#define _ewald_util_h

static char *SRCID_ewald_util_h = "";

#include <math.h>
#include "typedefs.h"
#include "complex.h"


extern real ewald_LRcorrection(FILE *fp,t_nsborder *nsb,
			       t_commrec *cr,t_forcerec *fr,
			       real charge[],t_block *excl,rvec x[],
			       rvec box_size,matrix lrvir);
/* Calculate the self energy and forces
 * when using Ewald/PME long range electrostatics.
 * Part of this is a constant, it is computed only once and stored in
 * a local variable. The remainder is computed every step.
 * PBC is taken into account. (Erik L.) 
 */

extern real calc_ewaldcoeff(real rc,real dtol);
/* Determines the Ewald parameter, both for Ewald and PME */
#endif
