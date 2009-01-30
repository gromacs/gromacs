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

#ifndef _mdrun_h
#define _mdrun_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "network.h"
#include "tgroup.h"
#include "filenm.h"
#include "mshift.h"
#include "force.h"
#include "time.h"
#include "edsam.h"
#include "mdebin.h"
#include "vcm.h"
#include "vsite.h"
#include "pull.h"
#include "update.h"
#include "genborn.h"

#define MD_GLAS         (1<<1)
#define MD_POLARISE     (1<<2)
#define MD_IONIZE       (1<<3)
#define MD_RERUN        (1<<4)
#define MD_FFSCAN       (1<<6)
#define MD_SEPPOT       (1<<7)
#define MD_PARTDEC      (1<<9)
#define MD_DDBONDCHECK  (1<<10)
#define MD_DDBONDCOMM   (1<<11)
#define MD_CONFOUT      (1<<12)
#define MD_NOGSTAT      (1<<13)
#define MD_REPRODUCIBLE (1<<14)
#define MD_READ_RNG     (1<<15)
#define MD_APPENDFILES  (1<<16)
#define MD_READ_EKIN    (1<<17)
#define MD_STARTFROMCPT (1<<18)


enum {
  ddnoSEL, ddnoINTERLEAVE, ddnoPP_PME, ddnoCARTESIAN, ddnoNR
};

typedef time_t gmx_integrator_t(FILE *log,t_commrec *cr,
				int nfile,t_filenm fnm[],
				bool bVerbose,bool bCompact,
				gmx_vsite_t *vsite,gmx_constr_t constr,
				int stepout,
				t_inputrec *inputrec,
				gmx_mtop_t *mtop,t_fcdata *fcd,
				t_state *state,rvec f[],
				rvec buf[],t_mdatoms *mdatoms,
				t_nrnb *nrnb,gmx_wallcycle_t wcycle,
				gmx_edsam_t ed, 
				t_forcerec *fr,
				gmx_genborn_t *born,
				int repl_ex_nst,int repl_ex_seed,
				real cpt_period,real max_hours,
				unsigned long Flags,
				int *nsteps_done);

/* ROUTINES from md.c */

extern gmx_integrator_t do_md;

/* ROUTINES from minimize.c */

extern gmx_integrator_t do_steep;
/* Do steepest descents EM */

extern gmx_integrator_t do_cg;
/* Do conjugate gradient EM */

extern gmx_integrator_t do_lbfgs;
/* Do conjugate gradient L-BFGS */

extern gmx_integrator_t do_nm;
/* Do normal mode analysis */

extern gmx_integrator_t do_tpi;
/* Do test particle insertion */


/* ROUTINES from sim_util.c */
extern void do_pbc_first(FILE *log,matrix box,t_forcerec *fr,
			 t_graph *graph,rvec x[]);

extern void do_pbc_first_mtop(FILE *fplog,int ePBC,matrix box,
			      gmx_mtop_t *mtop,rvec x[]);

extern void do_pbc_mtop(FILE *fplog,int ePBC,matrix box,
			gmx_mtop_t *mtop,rvec x[]);

		     
/* ROUTINES from stat.c */		
extern void global_stat(FILE *log,
			t_commrec *cr,gmx_enerdata_t *enerd,
			tensor fvir,tensor svir,rvec mu_tot,
			t_inputrec *inputrec,
			gmx_ekindata_t *ekind,bool bSumEkinhOld,
			gmx_constr_t constr,t_vcm *vcm,
			int *nabnsb,real *chkpt,real *terminate);
/* Communicate statistics over cr->mpi_comm_mysim */

void write_traj(FILE *fplog,t_commrec *cr,
		int fp_trn,bool bX,bool bV,bool bF,
		int fp_xtc,bool bXTC,int xtc_prec,
		char *fn_cpt,bool bCPT,
		gmx_mtop_t *top_global,
		int eIntegrator,int simulation_part,int step,double t,
		t_state *state_local,t_state *state_global,
		rvec *f_local,rvec *f_global);
/* Routine that writes frames to trn, xtc and/or checkpoint.
 * Data is collected to the master node only when necessary.
 */

extern int do_per_step(int step,int nstep);
/* Return TRUE if io should be done */

extern int do_any_io(int step, t_inputrec *ir);

/* ROUTINES from sim_util.c */
 
extern void print_time(FILE *out,time_t start,int step,t_inputrec *ir);

extern time_t print_date_and_time(FILE *log,int pid,char *title);

extern void nstop_cm(FILE *log,t_commrec *cr,
		     int start,int nr_atoms,real mass[],rvec x[],rvec v[]);

extern void finish_run(FILE *log,t_commrec *cr,char *confout,
		       t_inputrec *inputrec,
		       t_nrnb nrnb[],gmx_wallcycle_t wcycle,
		       double nodetime,double realtime,int step,
		       bool bWriteStat);

extern void calc_dispcorr(FILE *fplog,t_inputrec *ir,t_forcerec *fr,int step,
			  int natoms,matrix box,real lambda,
			  tensor pres,tensor virial,real ener[]);
     

typedef enum
{
  LIST_SCALARS	=0001,
  LIST_INPUTREC	=0002,
  LIST_TOP	=0004,
  LIST_X	=0010,
  LIST_V	=0020,
  LIST_F	=0040,
  LIST_LOAD	=0100
} t_listitem;

extern void check_nnodes_top(char *fn,t_topology *top);
/* Reset the tpr file to work with one node if necessary */

extern void init_single(FILE *log,
                        t_inputrec *inputrec, char *tpbfile, gmx_mtop_t *mtop,
			t_state *state);
     /*
      * Allocates space for the topology (top), the coordinates x, the
      * velocities v, masses mass. Reads the parameters, topology,
      * coordinates and velocities from the file specified in tpbfile
      */

extern void init_parallel(FILE *log,char *tpxfile,t_commrec *cr,
			  t_inputrec *inputrec,gmx_mtop_t *mtop,
			  t_state *state,
			  int list);
     /*
      * Loads the data for a simulation from the ring. Parameters, topology
      * coordinates, velocities, and masses are initialised equal to using
      * init_single() in the single processor version. The extra argument
      * f_add is allocated to use for the update of the forces, the load
      * array specifies in which part of the x and f array the subsystems
      * of the other processors are located. Homenr0, homenr1, nparts0 and
      * nparts1 are necessary to calculate the non bonded interaction using
      * the symmetry and thus calculating every force only once. List is a facility
      * for logging (and debugging). One can decide to print none or a set of
      * selected parameters to the file specified by log. Parameters are
      * printed by or-ing the corresponding items from t_listitem. A 0 (zero)
      * specifies that nothing is to be printed on the file. The function
      * returns the number of shifts over the ring to perform to calculate
      * all interactions.
      */

extern void start_time(void);
/* Start timing routines */

extern void update_time(void);
/* Update the timer.This must be done at least every INT_MAX microseconds,
 * or 2400 s, in order to give reliable answers.
 */
 
extern double node_time(void);
/* Return the node time so far in seconds. */

extern void do_shakefirst(FILE *log,gmx_constr_t constr,
			  t_inputrec *inputrec,t_mdatoms *md,
			  t_state *state,rvec buf[],rvec f[],
			  t_graph *graph,t_commrec *cr,t_nrnb *nrnb,
			  t_forcerec *fr,t_idef *idef);
			  
extern void dynamic_load_balancing(bool bVerbose,t_commrec *cr,real capacity[],
				   int dimension,t_mdatoms *md,t_topology *top,
				   rvec x[],rvec v[],matrix box);
/* Perform load balancing, i.e. split the particles over processors
 * based on their coordinates in the "dimension" direction.
 */
				   
extern void mdrunner(FILE *fplog,t_commrec *cr,int nfile,t_filenm fnm[],
		     bool bVerbose,bool bCompact,
		     ivec ddxyz,int dd_node_order,real rdd,real rconstr,
		     char *dddlb_opt,real dlb_scale,
		     char *ddcsx,char *ddcsy,char *ddcsz,
		     int nstepout,
		     gmx_edsam_t ed,int repl_ex_nst,int repl_ex_seed,
		     real pforce,real cpt_period,real max_hours,
		     unsigned long Flags);
/* Driver routine, that calls the different methods */

extern void init_md(FILE *fplog,
		    t_commrec *cr,t_inputrec *ir,
		    double *t,double *t0,
		    real *lambda,double *lam0,
		    t_nrnb *nrnb,gmx_mtop_t *mtop,
		    gmx_stochd_t *stochd,
		    int nfile,t_filenm fnm[],
		    int *fp_trn,int *fp_xtc,int *fp_ene,char **fn_cpt,
		    FILE **fp_dgdl,FILE **fp_field,
		    t_mdebin **mdebin,
		    tensor force_vir,tensor shake_vir,
		    rvec mu_tot,
		    bool *bNEMD,bool *bSimAnn,t_vcm **vcm, unsigned long Flags);
/* Routine in sim_util.c */
     
#endif	/* _mdrun_h */
