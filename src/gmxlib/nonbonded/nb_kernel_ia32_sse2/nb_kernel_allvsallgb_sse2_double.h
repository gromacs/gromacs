
#ifndef _NB_KERNEL_ALLVSALLGB_SSE2_DOUBLE_H
#define _NB_KERNEL_ALLVSALLGB_SSE2_DOUBLE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "types/simple.h"
#include "typedefs.h"

void
nb_kernel_allvsallgb_sse2_double(t_forcerec *           fr,
                                 t_mdatoms *            mdatoms,
                                 t_blocka *             excl,    
                                 double *               x,
                                 double *               f,
                                 double *               Vc,
                                 double *               Vvdw,
                                 double *               Vpol,
                                 int *                  outeriter,
                                 int *                  inneriter,
                                 void *                 work);

#endif
