
#ifndef _NB_KERNEL_ALLVSALL_H
#define _NB_KERNEL_ALLVSALL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "types/simple.h"
#include "typedefs.h"
#include "../nb_kernel.h"

void
nb_kernel_allvsall(t_nblist         * nlist,
                   rvec             * x,
                   rvec             * f,
                   t_forcerec       * fr,
                   t_mdatoms        * mdatoms,
                   nb_kernel_data_t * kernel_data,
                   t_nrnb           * nrnb);


#endif
