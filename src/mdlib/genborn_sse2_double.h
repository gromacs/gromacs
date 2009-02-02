#ifndef _genborn_sse2_double_h
#define _genborn_sse2_double_h

#include "typedefs.h"
#include "grompp.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

real 
calc_gb_chainrule_sse2_double(int natoms, t_nblist *nl, real *dadx, real *dvda, real *xd, real *f, int gb_algorithm, 
							  gmx_genborn_t *born);	

int 
calc_gb_rad_still_sse2_double(t_commrec *cr, t_forcerec *fr,int natoms, gmx_mtop_t *mtop,
							  const t_atomtypes *atype, real *x, t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md);

int 
calc_gb_rad_hct_sse2_double(t_commrec *cr, t_forcerec *fr, int natoms, gmx_mtop_t *mtop, 
							const t_atomtypes *atype, real *x, t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md);

int 
calc_gb_rad_obc_sse2_double(t_commrec *cr, t_forcerec * fr, int natoms, gmx_mtop_t *mtop,
							const t_atomtypes *atype, real *x, t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md);


#endif /* _genborn_sse2_double_h */