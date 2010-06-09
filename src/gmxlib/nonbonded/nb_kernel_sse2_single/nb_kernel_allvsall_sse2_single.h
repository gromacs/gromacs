
#ifndef _NB_KERNEL_ALLVSALL_SSE2_SINGLE_H
#define _NB_KERNEL_ALLVSALL_SSE2_SINGLE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "types/simple.h"
#include "typedefs.h"

void
nb_kernel_allvsall_sse2_single(t_forcerec *           fr,
                               t_mdatoms *            mdatoms,
                               t_blocka *             excl,    
                               real *                 x,
                               real *                 f,
                               real *                 Vc,
                               real *                 Vvdw,
                               int *                  outeriter,
                               int *                  inneriter,
                               void *                 work);

#endif
