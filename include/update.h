/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _update_h
#define _update_h

static char *SRCID_update_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) update.h 1.33 24 Jun 1996"
#endif /* HAVE_IDENT */

#include "typedefs.h"
#include "mshift.h"
#include "tgroup.h"
#include "network.h"
#include "force.h"
#include "pull.h"

extern void update(int          natoms,	/* number of atoms in simulation */
		   int      	start,
		   int          homenr,	/* number of home particles 	*/
		   int          step,
		   real         lambda, /* FEP scaling parameter */
		   real         *dvdlambda, /* FEP stuff */
		   t_parm       *parm,    /* input record and box stuff	*/
		   real         SAfactor, /* simulated annealing factor   */
		   t_mdatoms    *md,
		   rvec         x[],	/* coordinates of home particles */
		   t_graph      *graph,	
		   rvec         force[],/* forces on home particles 	*/
		   rvec         delta_f[],
		   rvec         vold[],	/* Old velocities		   */
		   rvec         vt[], 	/* velocities at whole steps */
		   rvec         v[],  	/* velocity at next halfstep   	*/
		   t_topology   *top,
		   t_groups     *grps,
		   tensor       vir_part,
		   t_commrec    *cr,
		   t_nrnb       *nrnb,
		   bool         bTYZ,
		   bool         bDoUpdate,
		   t_edsamyn    *edyn,
		   t_pull       *pulldata,
		   bool         bNEMD);
/* Return TRUE if OK, FALSE in case of Shake Error */
     
extern void calc_ke_part(bool bFirstStep,bool bSD,int start,int homenr,
			 rvec vold[],rvec v[],rvec vt[],
			 t_grpopts *opts,t_mdatoms *md,
			 t_groups *grps,t_nrnb *nrnb,
			 real lambda,real *dvdlambda);
/*
 * Compute the partial kinetic energy for home particles;
 * will be accumulated in the calling routine.
 * The velocity a the whole step is obtained by averaging
 * the velocities of minus and plus a half step,
 * in case of Stochastic Dynamics a correction is applied.
 * The tensor is
 *
 * Ekin = SUM(i) 0.5 m[i] v[i] (x) v[i]
 *    
 *     use v[i] = v[i] - u[i] when calculating temperature
 *
 * u must be accumulated already.
 *
 * Now also computes the contribution of the kinetic energy to the
 * free energy
 *
 */

extern void calc_ke_part_visc(bool bFirstStep,int start,int homenr,
			      matrix box,rvec x[],
			      rvec vold[],rvec v[],rvec vt[],
			      t_grpopts *opts,t_mdatoms *md,
			      t_groups *grps,t_nrnb *nrnb,
			      real lambda,real *dvdlambda);
/* The same as calc_ke_part, but for viscosity calculations.
 * The cosine velocity profile is excluded from the kinetic energy.
 * The new amplitude of the velocity profile is calculated for this
 * node and stored in grps->cosacc.mvcos.
 */

extern void init_sd_consts(int ngtc,real tau_t[],real dt);
/* Initialization of the SD constants (obviously). */

/* Routines from coupling.c to do with Temperature, Pressure and coupling
 * algorithms.
 */
extern real run_aver(real old,real cur,int step,int nmem);

extern void berendsen_tcoupl(t_grpopts *opts,t_groups *grps,
			     real dt,real SAfactor);
extern void nosehoover_tcoupl(t_grpopts *opts,t_groups *grps,
			      real dt,real SAfactor);
/* Compute temperature scaling. For Nose-Hoover it is done in update. */

extern real calc_temp(real ekin,real nrdf);
/* Calculate the temperature */

extern void calc_pres(int ePBC,matrix box,
		      tensor ekin,tensor vir,tensor pres,real Elr);
/* Calculate the pressure. Unit of pressure is bar, If Elr != 0
 * a long range correction based on Ewald/PPPM is made (see c-code)
 */

extern void parrinellorahman_pcoupl(t_inputrec *ir,int step,tensor pres,
				   tensor box,tensor boxv,tensor M);
  
extern void berendsen_pcoupl(t_inputrec *ir,int step,tensor pres,
			     matrix box,int start,int nr_atoms,
			     rvec x[],unsigned short cFREEZE[],
			     t_nrnb *nrnb,ivec nFreeze[]);

extern void correct_ekin(FILE *log,int start,int end,rvec v[],
			 rvec vcm,real mass[],real tmass,tensor ekin);
/* Correct ekin for vcm */

#endif	/* _update_h */

