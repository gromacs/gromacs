/*
 * $Id$
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

#ifndef _update_h
#define _update_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "mshift.h"
#include "tgroup.h"
#include "network.h"
#include "force.h"
#include "pull.h"
#include "gmx_random.h"

/* Abstract type for stochastic dynamics */
typedef struct gmx_stochd *gmx_stochd_t;

/* Initialize the stochastic dynamics struct */
extern gmx_stochd_t init_stochd(FILE *fplog,t_inputrec *ir);

/* Store the random state from sd in state */
extern void get_stochd_state(gmx_stochd_t sd,t_state *state);

/* Set the random in sd from state */
extern void set_stochd_state(gmx_stochd_t sd,t_state *state);

/* Store the box at step step
 * as a reference state for simulations with box deformation.
 */
extern void set_deform_reference_box(int step,matrix box);

extern void update(FILE         *fplog,
		   int          step,
		   real         *dvdlambda, /* FEP stuff */
		   t_inputrec   *inputrec,  /* input record and box stuff	*/
		   t_mdatoms    *md,
		   t_state      *state,
		   t_graph      *graph,	
		   rvec         force[],    /* forces on home particles */
		   rvec         xprime[],   /* buffer for x for update  */
		   t_fcdata     *fcd,
		   t_topology   *top,
		   t_groups     *grps,
		   tensor       vir_part,
		   matrix       *scale_tot,
		   t_commrec    *cr,
		   t_nrnb       *nrnb,
		   gmx_wallcycle_t wcycle,
		   gmx_stochd_t sd,
		   gmx_constr_t constr,
		   bool         bHaveConstr,
		   bool         bNEMD,
		   bool         bInitStep);
/* Return TRUE if OK, FALSE in case of Shake Error */
     
extern void calc_ke_part(rvec v[],t_grpopts *opts,t_mdatoms *md,
			 t_groups *grps,t_nrnb *nrnb,
			 real lambda);
/*
 * Compute the partial kinetic energy for home particles;
 * will be accumulated in the calling routine.
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

extern void calc_ke_part_visc(matrix box,rvec x[],rvec v[],
			      t_grpopts *opts,t_mdatoms *md,
			      t_groups *grps,t_nrnb *nrnb,
			      real lambda);
/* The same as calc_ke_part, but for viscosity calculations.
 * The cosine velocity profile is excluded from the kinetic energy.
 * The new amplitude of the velocity profile is calculated for this
 * node and stored in grps->cosacc.mvcos.
 */

extern real run_aver(real old,real cur,int step,int nmem);

extern void berendsen_tcoupl(t_grpopts *opts,t_groups *grps,real dt);

extern void nosehoover_tcoupl(t_grpopts *opts,t_groups *grps,real dt,
			      real xi[],double ixi[]);
/* Compute temperature scaling. For Nose-Hoover it is done in update. */

extern real nosehoover_energy(t_grpopts *opts,t_groups *grps,
			      real *xi,double *ixi);
/* Returns the Nose-Hoover contribution to the conserved energy */

extern void vrescale_tcoupl(t_grpopts *opts,t_groups *grps,real dt,
			    double therm_integral[],
			    gmx_rng_t rng);
/* Compute temperature scaling. For V-rescale it is done in update. */

extern real vrescale_energy(t_grpopts *opts,double therm_integral[]);
/* Returns the V-rescale contribution to the conserved energy */

/* Set reference temp for simulated annealing at time t*/
extern void update_annealing_target_temp(t_grpopts *opts,real t); 

extern real calc_temp(real ekin,real nrdf);
/* Calculate the temperature */

extern real calc_pres(int ePBC,int nwall,matrix box,
		      tensor ekin,tensor vir,tensor pres,real Elr);
/* Calculate the pressure tensor, returns the scalar pressure.
 * The unit of pressure is bar, If Elr != 0
 * a long range correction based on Ewald/PPPM is made (see c-code)
 */

extern void parrinellorahman_pcoupl(FILE *fplog,int step,
				    t_inputrec *ir,tensor pres,
				    tensor box,tensor box_rel,tensor boxv,
				    tensor M,matrix *scale_tot,
				    bool bFirstStep);
  
extern void berendsen_pcoupl(FILE *fplog,int step,
			     t_inputrec *ir,tensor pres,matrix box,
			     matrix mu);


extern void berendsen_pscale(t_inputrec *ir,matrix mu,
			     matrix box,matrix box_rel,
			     int start,int nr_atoms,
			     rvec x[],unsigned short cFREEZE[],
			     t_nrnb *nrnb);

extern void correct_ekin(FILE *log,int start,int end,rvec v[],
			 rvec vcm,real mass[],real tmass,tensor ekin);
/* Correct ekin for vcm */

#endif	/* _update_h */

