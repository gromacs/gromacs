
#ifndef _NB_KERNEL_ALLVSALL_H
#define _NB_KERNEL_ALLVSALL_H

#include <types/simple.h>


		
void
nb_kernel_allvsall(t_forcerec *           fr,
				   t_mdatoms *            mdatoms,
				   t_blocka *             excl,    
				   real *                 x,
				   real *                 f,
				   real *                 Vc,
				   real *                 Vvdw,
				   int *                  outeriter,
				   int *                  inneriter,
				   void *                 work);

void
nb_kernel_allvsall_gb(t_forcerec *           fr,
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