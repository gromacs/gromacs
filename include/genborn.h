/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */


#ifndef _genborn_h
#define _genborn_h

#include "typedefs.h"
#include "grompp.h"


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


typedef struct
{
	int nbonds;
	int bond[10];
	real length[10];
} genborn_bonds_t;


/* Struct to hold all the information for GB 
 * All these things are currently allocated in md.c
 */
typedef struct
{
	int nr;                   /* number of atoms, length of arrays below */
	int n12;                  /* number of 1-2 (bond) interactions       */
	int n13;                  /* number of 1-3 (angle) terms             */
	int n14;                  /* number of 1-4 (torsion) terms           */
	int nlocal;               /* Length of local arrays (with DD)        */
	
	/* Arrays below that end with _globalindex are used for setting up initial values of
	 * all gb parameters and values. They all have length natoms, which for DD is the 
	 * global atom number. 
	 * Values are then taken from these arrays to local copies, that have names without
	 * _globalindex, in the routine make_local_gb(), which is called once for single
	 * node runs, and for DD at every call to dd_partition_system
	 */

	real  *gpol;              /* Atomic polarisation energies */
	real  *gpol_globalindex;  /*  */
	real  *gpol_still_work;   /* Work array for Still model */
	real  *gpol_hct_work;     /* Work array for HCT/OBC models */
	real  *bRad;              /* Atomic Born radii */
	real  *vsolv;             /* Atomic solvation volumes */
	real  *vsolv_globalindex; /*  */
	
	int *vs;                  /* Array for vsites-exclusions */   
	int *vs_globalindex;      /*  */
		
	real es;                  /* Solvation energy and derivatives */
	real *asurf;              /* Atomic surface area */
	rvec *dasurf;             /* Surface area derivatives */
	real as;                  /* Total surface area */

	real *drobc;              /* Parameters for OBC chain rule calculation */
	real *param;              /* Precomputed factor rai*atype->S_hct for HCT/OBC */
	real *param_globalindex;  /*  */
	
	real *log_table;          /* Table for logarithm lookup */
	
	real obc_alpha;           /* OBC parameters */
	real obc_beta;            /* OBC parameters */
	real obc_gamma;           /* OBC parameters */
	real gb_doffset;          /* Dielectric offset for Still/HCT/OBC */
	
	real *work;               /* Used for parallel summation and in the chain rule, length natoms         */
	real *dd_work;            /* Used for domain decomposition parallell runs, length natoms              */
	int  *count;              /* Used for setting up the special gb nblist, length natoms                 */
	int  **nblist_work;       /* Used for setting up the special gb nblist, dim natoms*nblist_work_nalloc */
	int  nblist_work_nalloc;  /* Length of second dimension of nblist_work                                */
} 
gmx_genborn_t;


/* Initialise GB stuff */
int init_gb(gmx_genborn_t **p_born,t_commrec *cr, t_forcerec *fr, t_inputrec *ir,
			gmx_mtop_t *mtop, rvec x[], real rgbradii, int gb_algorithm);


/* Born radii calculations, both with and without SSE acceleration */
int calc_gb_rad(t_commrec *cr, t_forcerec *fr, t_inputrec *ir,gmx_localtop_t *top,
				const t_atomtypes *atype, rvec x[], rvec f[],t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md);



/* Bonded GB interactions */								
real gb_bonds_tab(real *x, real *f, real *charge, real *p_gbtabscale,
				  real *invsqrta, real *dvda, real *GBtab, t_idef *idef,
				  real epsilon_r, real facel);



void gb_pd_send(t_commrec *cr, real *send_data, int nr);


/* Functions for setting up the F_GB12,13,14 lists in grompp */
int 
init_gb_plist(t_params *p_list);

int 
convert_gb_params(gmx_ffparams_t *ffparams, t_functype ftype, t_params *gb_plist, t_ilist *il);

int 
generate_gb_topology(gmx_mtop_t *mtop, t_molinfo *mi);
						 



/* Functions for calculating adjustments due to ie chain rule terms */
real 
calc_gb_forces(t_commrec *cr, t_mdatoms *md, gmx_genborn_t *born, gmx_localtop_t *top, const t_atomtypes *atype,
			   rvec x[], rvec f[], t_forcerec *fr,t_idef *idef,int gb_algorithm, bool bRad);


int
make_gb_nblist(t_commrec *cr, int natoms, int gb_algorithm, real gbcut, rvec x[], 
			   t_forcerec *fr, t_idef *idef, gmx_genborn_t *born);

void 
make_local_gb(t_commrec *cr, gmx_genborn_t *born, int gb_algorithm);


#endif /* _genborn_h */

