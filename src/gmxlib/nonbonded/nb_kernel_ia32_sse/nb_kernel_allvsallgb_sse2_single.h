
#ifndef _NB_KERNEL_ALLVSALLGB_SSE2_SINGLE_H
#define _NB_KERNEL_ALLVSALLGB_SSE2_SINGLE_H

#include "types/simple.h"
#include "typedefs.h"

void
nb_kernel_allvsallgb_sse2_single(t_forcerec *           fr,
                                 t_mdatoms *            mdatoms,
                                 t_blocka *             excl,    
                                 float *                x,
                                 float *                f,
                                 float *                Vc,
                                 float *                Vvdw,
                                 float *                Vpol,
                                 int *                  outeriter,
                                 int *                  inneriter,
                                 void *                 work);

#endif
