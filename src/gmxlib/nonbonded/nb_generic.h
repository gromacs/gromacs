#ifndef _nb_generic_h_
#define _nb_generic_h_

#include "types/simple.h"
#include "typedefs.h"

void
gmx_nb_generic_kernel(t_nblist *           nlist,
					  t_forcerec *         fr,
					  t_mdatoms *          mdatoms,
					  real *               x,
					  real *               f,
					  real *               fshift,
					  real *               Vc,
					  real *               Vvdw,
					  real                 tabscale,  
					  real *               VFtab,
					  int *                outeriter,
					  int *                inneriter);

#endif

