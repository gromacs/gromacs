/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_inloop_h = "$Id$";

 
#ifndef _in_loop_h_
#define _in_loop_h_

#include "typedefs.h"
	
/******************************************************************
 *
 *        C    I N N E R   L O O P S
 *
 ******************************************************************/
extern void c_coul(real ix,real iy,real iz,real qi,
		   real pos[],int nj,t_nl_j jjnr[],
		   real charge[],
		   real faction[],real fip[3],
		   real *Vc);
		   
extern void c_ljc(real ix,real iy,real iz,real qi,
		  real pos[],int nj,int type[],t_nl_j jjnr[],
		  real charge[],real nbfp[],
		  real faction[],real fip[],
		  real *Vc,real *Vnb);
		  
extern void c_bham(real ix,real iy,real iz,real qi,
		   real pos[],int nj,int type[],t_nl_j jjnr[],
		   real charge[],real nbfp[],
		   real faction[],real fip[],
		   real *Vc,real *Vnb);

extern void c_water(int  i0,real xw[],real eps,
		    real pos[],int nj,int type[],t_nl_j jjnr[],
		    real charge[],real nbfp[],
		    real faction[],real fw[],
		    real *Vc,real *Vnb);
/* Water loop */

extern void c_wcoul(int  i0,real xw[],real eps,
		    real pos[],int nj,t_nl_j jjnr[],
		    real charge[],real faction[],real fw[],real *Vc);
/* Water loop including only coulomb */

extern void c_ljcfree(real ix,real iy,real iz,int inr,
		      real pos[],int nj,t_nl_j jjnr[],
		      int  typeA[],  int typeB[],
		      real eps,
		      real chargeA[],real chargeB[],
		      real nbfpA[],  real nbfpB[],
		      real faction[],real fip[],
		      real *Vc,      real *Vnb,
		      real lambda,   real *dvdlambda,
		      real k_rf,     real c_rf,
		      real tfac,real trunctab[]);
/* LJ+Coulomb with free energy perturbation calculation 
 * This also contains a reaction field term. If rffac is zero,
 * no RF is applied.
 */

extern void c_free(real ix,real iy,real iz,int inr,
		   real pos[],int nj,t_nl_j jjnr[],
		   int  typeA[],  int typeB[],
		   real epsilon,
		   real chargeA[],real chargeB[],
		   real nbfpA[],  real nbfpB[],
		   real faction[],real fip[],
		   real *Vc,      real *Vnb,
		   real lambda,   real *dvdlambda,
		   int  ntab,     real tabscale,
		   real VFtab[]);
/* LJ+Coulomb with free energy perturbation calculation 
 * Potential and forces are calculated using tables.
 */

extern void c_tab(real ix,real iy,real iz,real qi,
		  real pos[],int nj,int type[],int jjnr[],real charge[],
		  real nbfp[],real faction[],real fip[],
		  real *Vc,real *Vnb,int ntab,real tabscale,real VFtab[]);

extern void c_bhamtab(real ix,real iy,real iz,real qi,
		      real pos[],int nj,int type[],int jjnr[],real charge[],
		      real nbfp[],real faction[],real fip[],
		      real *Vc,real *Vnb,int ntab,
		      real tabscale,real tabscale_bham,real VFtab[]);

extern void c_coultab(real ix,real iy,real iz,real qi,
		      real pos[],int nj,int jjnr[],real charge[],
		      real faction[],real fip[],
		      real *Vc,int ntab,real tabscale,real VFtab[]);
/* These two routines calculate force based on tables.
 * The first does both LJ and Coulomn, while the second only does coulomb.
 * More info in inloopc.c
 */
extern void c_watertab(int i0,real xw[],real fudge,
		       real pos[],int nj,int type[],t_nl_j jjnr[],
		       real charge[],real nbfp[],
		       real faction[],real fw[],
		       real *Vc,real *Vnb,int ntab,real tabscale,real VFtab[]);
/* Water loop including tables */

void c_wcoultab(int i0,real xw[],real fudge,
		real pos[],int nj,t_nl_j jjnr[],
		real charge[],real faction[],real fw[],
		real *Vc,int ntab,real tabscale,real VFtab[]);
/*  Water loop including only coulomb, but with a table function */

#endif /* _in_loop_h_ */
