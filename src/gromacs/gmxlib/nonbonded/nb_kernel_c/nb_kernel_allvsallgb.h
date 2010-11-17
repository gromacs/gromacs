
#ifndef _NB_KERNEL_ALLVSALLGB_H
#define _NB_KERNEL_ALLVSALLGB_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "types/simple.h"
#include "typedefs.h"

void
nb_kernel_allvsallgb(t_forcerec *           fr,
                     t_mdatoms *            mdatoms,
                     t_blocka *             excl,    
                     real *                 x,
                     real *                 f,
                     real *                 Vc,
                     real *                 Vvdw,
                     real *                 vpol,
                     int *                  outeriter,
                     int *                  inneriter,
                     void *                 work);

#endif
