/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _bondf_h
#define _bondf_h


#include <stdio.h>
#include "typedefs.h"
#include "nrnb.h"
#include "pbc.h"
#include "genborn.h"

#ifdef __cplusplus
extern "C" {
#endif

int glatnr(int *global_atom_index,int i);
/* Returns the global topology atom number belonging to local atom index i.
 * This function is intended for writing ascii output
 * and returns atom numbers starting at 1.
 * When global_atom_index=NULL returns i+1.
 */

void calc_bonds(FILE *fplog,const gmx_multisim_t *ms,
                const t_idef *idef,
                rvec x[],history_t *hist,
                rvec f[],t_forcerec *fr,
                const t_pbc *pbc,const t_graph *g,
                gmx_enerdata_t *enerd,t_nrnb *nrnb,real lambda,
                const t_mdatoms *md,
                t_fcdata *fcd,int *ddgatindex,
                t_atomtypes *atype, gmx_genborn_t *born,
                gmx_bool bPrintSepPot,gmx_large_int_t step);
/* 
 * The function calc_bonds() calculates all bonded force interactions.
 * The "bonds" are specified as follows:
 *   int nbonds
 *	    the total number of bonded interactions.
 *   t_iatom *forceatoms
 *     specifies which atoms are involved in a bond of a certain 
 *     type, see also struct t_idef.
 *   t_functype *functype
 *	    defines for every bonded force type what type of function to 
 *     use, see also struct t_idef.
 *   t_iparams *forceparams
 *	    defines the parameters for every bond type, see also struct 
 *     t_idef.
 *   real epot[NR_F]
 *     total potential energy split up over the function types.
 *   int *ddgatindex
 *     global atom number indices, should be NULL when not using DD.
 *   gmx_bool bPrintSepPot
 *     if TRUE print local potential and dVdlambda for each bonded type.
 *   int step
 *     used with bPrintSepPot
 *   return value:
 *	    the total potential energy (sum over epot).
 */

void calc_bonds_lambda(FILE *fplog,
			      const t_idef *idef,
			      rvec x[],
			      t_forcerec *fr,
			      const t_pbc *pbc,const t_graph *g,
			      gmx_enerdata_t *enerd,t_nrnb *nrnb,
			      real lambda,
			      const t_mdatoms *md,
			      t_fcdata *fcd,int *global_atom_index);
/* As calc_bonds, but only determines the potential energy
 * for the perturbed interactions.
 * The shift forces in fr are not affected.
 */

real posres(int nbonds,
		   const t_iatom forceatoms[],const t_iparams forceparams[],
		   const rvec x[],rvec f[],rvec vir_diag,
		   t_pbc *pbc,
		   real lambda,real *dvdlambda,
		   int refcoord_scaling,int ePBC,rvec comA,rvec comB);
/* Position restraints require a different pbc treatment from other bondeds */

real bond_angle(const rvec xi,const rvec xj,const rvec xk,
		       const t_pbc *pbc,
		       rvec r_ij,rvec r_kj,real *costh,
		       int *t1,int *t2);	/* out */
/* Calculate bond-angle. No PBC is taken into account (use mol-shift) */

real dih_angle(const rvec xi,const rvec xj,const rvec xk,const rvec xl,
		      const t_pbc *pbc,
		      rvec r_ij,rvec r_kj,rvec r_kl,rvec m,rvec n, /* out */
		      real *sign,
		      int *t1,int *t2,int *t3);
/* Calculate dihedral-angle. No PBC is taken into account (use mol-shift) */

void do_dih_fup(int i,int j,int k,int l,real ddphi,
		       rvec r_ij,rvec r_kj,rvec r_kl,
		       rvec m,rvec n,rvec f[],rvec fshift[],
		       const t_pbc *pbc,const t_graph *g,
		       const rvec *x,int t1,int t2,int t3);
/* Do an update of the forces for dihedral potentials */

/*************************************************************************
 *
 *  Bonded force functions
 *
 *************************************************************************/
  t_ifunc bonds,g96bonds,morse_bonds,cubic_bonds,FENE_bonds,restraint_bonds;
  t_ifunc angles,g96angles,cross_bond_bond,cross_bond_angle,urey_bradley,quartic_angles;
  t_ifunc pdihs,idihs,rbdihs;
  t_ifunc tab_bonds,tab_angles,tab_dihs;
  t_ifunc polarize,water_pol,thole_pol,angres,angresz,unimplemented;

#ifdef __cplusplus
}
#endif

#endif	/* _bondf_h */
