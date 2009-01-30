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
} bonds_t;

typedef struct
{
  real length[20000];
  real angle[20000];
} bl_t;

/* Struct to hold all the information for GB 
 * All these things are currently allocated in md.c
 */
typedef struct
{
	int nr;                  /* number of atoms, length of arrays below */
	int n12;                 /* number of 1-2 (bond) interactions       */
	int n13;                 /* number of 1-3 (angle) terms             */
	int n14;                 /* number of 1-4 (torsion) terms           */
	real  *gpol;             /* Atomic polarisation energies */
	real  *bRad;             /* Atomic Born radii */
	real  *vsolv;            /* Atomic solvation volumes */
	int *vs;                 /* Array for vsites-exclusions */     
	int nvs;                 /* Length of vs array         */
	real es;                 /* Solvation energy and derivatives */
	real *asurf;             /* Atomic surface area */
	rvec *dasurf;            /* Surface area derivatives */
	real as;                 /* Total surface area */
	real *S_hct;             /* Overlap factors for HCT method */
	real *drobc;             /* Parameters for OBC chain rule calculation */
	real *param;             /* Precomputed factor raj*atype->S_hct for HCT/OBC */
	real *log_table;         /* Table for logarithm lookup */
	
	real obc_alpha;          /* OBC parameters */
	real obc_beta;           /* OBC parameters */
	real obc_gamma;          /* OBC parameters */
	real gb_doffset;         /* Dielectric offset for Still/HCT/OBC */
	
	real *work;              /* used for parallel summation */
} 
gmx_genborn_t;


/* Initialise GB stuff */
int init_gb(gmx_genborn_t **p_born,t_commrec *cr, t_forcerec *fr, t_inputrec *ir,
			gmx_mtop_t *mtop, rvec x[], real rgbradii, int gb_algorithm);


/* Born radii calculations, both with and without SSE acceleration */
int calc_gb_rad(t_commrec *cr, t_forcerec *fr, t_inputrec *ir,int natoms, int nrfa, gmx_mtop_t *mtop,
				const t_atomtypes *atype, rvec x[], rvec f[],t_nblist *nl, gmx_genborn_t *born,t_mdatoms *md);



/* Bonded GB interactions */								
real gb_bonds_tab(int nbonds, real *x, real *f, real *charge, real *p_gbtabscale,
				  real *invsqrta, real *dvda, real *GBtab, const t_iatom forceatoms[],
				  real epsilon_r, real facel);



void gb_pd_send(t_commrec *cr, real *send_data, int nr);


/* Functions for setting up the F_GB list in grompp */
int 
init_gb_plist(t_params *p_list);

int 
convert_gb_params(t_idef *idef, t_functype ftype, int start, t_params *gb_plist, gmx_genborn_t *born);

int 
generate_gb_topology(gmx_mtop_t *mtop, t_params *plist, t_params *gb_plist, gmx_genborn_t *born);
						 



/* Functions for calculating adjustments due to ie chain rule terms */
real 
calc_gb_forces(t_commrec *cr, t_mdatoms *md, gmx_genborn_t *born, gmx_mtop_t *mtop, const t_atomtypes *atype, int nr, 
			   rvec x[], rvec f[], t_forcerec *fr,const t_iatom forceatoms[],int gb_algorithm, bool bRad);


int
gb_nblist_siev(t_commrec *cr, int natoms, int gb_algorithm, real gbcut, rvec x[], t_forcerec *fr, t_ilist *il, int n14);


#endif /* _genborn_h */

