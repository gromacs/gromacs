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

#ifndef _update_h
#define _update_h

#include "typedefs.h"
#include "mshift.h"
#include "tgroup.h"
#include "network.h"
#include "force.h"
#include "pull.h"
#include "gmx_random.h"
#include "maths.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Abstract type for stochastic dynamics */
typedef struct gmx_update *gmx_update_t;

/* Initialize the stochastic dynamics struct */
gmx_update_t init_update(FILE *fplog,t_inputrec *ir);

/* Store the random state from sd in state */
void get_stochd_state(gmx_update_t sd,t_state *state);

/* Set the random in sd from state */
void set_stochd_state(gmx_update_t sd,t_state *state);

/* Store the box at step step
 * as a reference state for simulations with box deformation.
 */
void set_deform_reference_box(gmx_update_t upd,
				     gmx_large_int_t step,matrix box);

void update_tcouple(FILE         *fplog,
			   gmx_large_int_t   step,
			   t_inputrec   *inputrec,   
			   t_state      *state,
			   gmx_ekindata_t *ekind,
			   gmx_wallcycle_t wcycle,
			   gmx_update_t upd,
			   t_extmass    *MassQ,
			   t_mdatoms    *md
);

void update_pcouple(FILE         *fplog,
			   gmx_large_int_t   step,
			   t_inputrec   *inputrec,   
			   t_state      *state,
			   matrix       pcoupl_mu,
			   matrix       M,
			   gmx_wallcycle_t wcycle,
			   gmx_update_t upd,
			   gmx_bool         bInitStep);

void update_coords(FILE         *fplog,
			  gmx_large_int_t   step,
			  t_inputrec   *inputrec,  /* input record and box stuff	*/
			  t_mdatoms    *md,
			  t_state      *state,
			  rvec         *f,    /* forces on home particles */
			  gmx_bool         bDoLR,
			  rvec         *f_lr,
			  t_fcdata     *fcd,
			  gmx_ekindata_t *ekind,
			  matrix       M,
			  gmx_wallcycle_t wcycle,
			  gmx_update_t upd,
			  gmx_bool         bInitStep,
			  int          bUpdatePart,
			  t_commrec    *cr,  /* these shouldn't be here -- need to think about it */
			  t_nrnb       *nrnb,
			  gmx_constr_t constr,
			  t_idef       *idef);

/* Return TRUE if OK, FALSE in case of Shake Error */

void update_constraints(FILE         *fplog,
			       gmx_large_int_t   step,
			       real         *dvdlambda, /* FEP stuff */
			       t_inputrec   *inputrec,  /* input record and box stuff	*/
			       gmx_ekindata_t *ekind,
			       t_mdatoms    *md,
			       t_state      *state,
			       t_graph      *graph,	
			       rvec         force[],    /* forces on home particles */
			       t_idef       *idef,
			       tensor       vir_part,
			       tensor       vir,
			       t_commrec    *cr,
			       t_nrnb       *nrnb,
			       gmx_wallcycle_t wcycle,
			       gmx_update_t upd,
			       gmx_constr_t constr,
			       gmx_bool         bInitStep,
			       gmx_bool         bFirstHalf,
			       gmx_bool         bCalcVir,
			       real         vetanew);

/* Return TRUE if OK, FALSE in case of Shake Error */

void update_box(FILE         *fplog,
		       gmx_large_int_t   step,
		       t_inputrec   *inputrec,  /* input record and box stuff	*/
		       t_mdatoms    *md,
		       t_state      *state,
		       t_graph      *graph,
		       rvec         force[],    /* forces on home particles */
		       matrix       *scale_tot,
		       matrix       pcoupl_mu,		    
		       t_nrnb       *nrnb,
		       gmx_wallcycle_t wcycle,
		       gmx_update_t upd,
		       gmx_bool         bInitStep,
		       gmx_bool         bFirstHalf);
/* Return TRUE if OK, FALSE in case of Shake Error */

void calc_ke_part(t_state *state,t_grpopts *opts,t_mdatoms *md,
			 gmx_ekindata_t *ekind,t_nrnb *nrnb,gmx_bool bEkinAveVel, gmx_bool bSaveOld);
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

void
init_ekinstate(ekinstate_t *ekinstate,const t_inputrec *ir);

void
update_ekinstate(ekinstate_t *ekinstate,gmx_ekindata_t *ekind);

void
restore_ekinstate_from_state(t_commrec *cr,
			     gmx_ekindata_t *ekind,ekinstate_t *ekinstate);

void berendsen_tcoupl(t_inputrec *ir,gmx_ekindata_t *ekind,real dt);

void nosehoover_tcoupl(t_grpopts *opts,gmx_ekindata_t *ekind,real dt,
			      double xi[],double vxi[],t_extmass *MassQ);

t_state *init_bufstate(const t_state *template_state);

void destroy_bufstate(t_state *state);

void trotter_update(t_inputrec *ir, gmx_large_int_t step, gmx_ekindata_t *ekind, 
			   gmx_enerdata_t *enerd, t_state *state, tensor vir, t_mdatoms *md, 
			   t_extmass *MassQ, int **trotter_seqlist, int trotter_seqno);

int **init_npt_vars(t_inputrec *ir, t_state *state, t_extmass *Mass, gmx_bool bTrotter); 

real NPT_energy(t_inputrec *ir, t_state *state, t_extmass *MassQ);
/* computes all the pressure/tempertature control energy terms to get a conserved energy */

void NBaroT_trotter(t_grpopts *opts, real dt,
			   double xi[],double vxi[],real *veta, t_extmass *MassQ);

void vrescale_tcoupl(t_inputrec *ir,gmx_ekindata_t *ekind,real dt,
			    double therm_integral[],
			    gmx_rng_t rng);
/* Compute temperature scaling. For V-rescale it is done in update. */

real vrescale_energy(t_grpopts *opts,double therm_integral[]);
/* Returns the V-rescale contribution to the conserved energy */

void rescale_velocities(gmx_ekindata_t *ekind,t_mdatoms *mdatoms,
			       int start,int end,rvec v[]);
/* Rescale the velocities with the scaling factor in ekind */

void update_annealing_target_temp(t_grpopts *opts,real t); 
/* Set reference temp for simulated annealing at time t*/

real calc_temp(real ekin,real nrdf);
/* Calculate the temperature */

real calc_pres(int ePBC,int nwall,matrix box,
		      tensor ekin,tensor vir,tensor pres,real Elr);
/* Calculate the pressure tensor, returns the scalar pressure.
 * The unit of pressure is bar, If Elr != 0
 * a long range correction based on Ewald/PPPM is made (see c-code)
 */

void parrinellorahman_pcoupl(FILE *fplog,gmx_large_int_t step,
				    t_inputrec *ir,real dt,tensor pres,
				    tensor box,tensor box_rel,tensor boxv,
				    tensor M,matrix mu,
				    gmx_bool bFirstStep);
  
void berendsen_pcoupl(FILE *fplog,gmx_large_int_t step,
			     t_inputrec *ir,real dt,tensor pres,matrix box,
			     matrix mu);


void berendsen_pscale(t_inputrec *ir,matrix mu,
			     matrix box,matrix box_rel,
			     int start,int nr_atoms,
			     rvec x[],unsigned short cFREEZE[],
			     t_nrnb *nrnb);

void correct_ekin(FILE *log,int start,int end,rvec v[],
			 rvec vcm,real mass[],real tmass,tensor ekin);
/* Correct ekin for vcm */


#ifdef __cplusplus
}
#endif

#endif	/* _update_h */

