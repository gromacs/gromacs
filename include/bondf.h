/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GROningen MAchine for Chemical Simulation
 */

#ifndef _bondf_h
#define _bondf_h

static char *SRCID_bondf_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) bondf.h 1.28 2/19/97"
#endif /* HAVE_IDENT */

#ifdef CPLUSPLUS
extern "C" {
#endif

#include <stdio.h>
#include "typedefs.h"
#include "nrnb.h"
#include "pbc.h"

extern void calc_bonds(FILE *log,t_idef *idef,
                       rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
                       real epot[],t_nrnb *nrnb,matrix box,real lambda,
		       t_mdatoms *md,int ngrp,real egnb[],real egcoul[]);
/* 
 * The function calc_bonds() caluclates all bonded force interactions.
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
 *   return value:
 *	    the total potential energy (sum over epot).
 */

extern real bond_angle(FILE *log,matrix box,
		       rvec xi,rvec xj,rvec xk,	                /* in  */
		       rvec r_ij,rvec r_kj,real *costh);	/* out */
/* Calculate bond-angle. No PBC is taken into account (use mol-shift) */

extern real dih_angle(FILE *log,matrix box,
		      rvec xi,rvec xj,rvec xk,rvec xl,   /* in */
		      rvec r_ij,rvec r_kj,rvec r_kl,rvec m,rvec n, /* out */
		      real *cos_phi,real *sign);
/* Calculate dihedral-angle. No PBC is taken into account (use mol-shift) */

extern void do_dih_fup(FILE *log,int i,int j,int k,int l,real ddphi,
		       rvec r_ij,rvec r_kj,rvec r_kl,
		       rvec m,rvec n,rvec f[],t_forcerec *fr,t_graph *g,
		       rvec x[]);
/* Do an update of the forces for dihedral potentials */

/*************************************************************************
 *
 *  Bonded force functions
 *
 *************************************************************************/
  extern t_ifunc bonds,g96bonds,morsebonds;
  extern t_ifunc angles,g96angles;
  extern t_ifunc pdihs,idihs,rbdihs;
  extern t_ifunc water_pol,posres,do_14,unimplemented;

#ifdef CPLUSPLUS
}
#endif

#endif	/* _bondf_h */
