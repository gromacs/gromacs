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
#include "nsb.h"
#include "mshift.h"
#include "force.h"
#include "time.h"
#include "edsam.h"
#include "mdebin.h"
#include "vcm.h"
#include "vsite.h"
#include "pull.h"

#define MD_MULTISIM  (1<<0)
#define MD_GLAS      (1<<1)
#define MD_POLARISE  (1<<2)
#define MD_IONIZE    (1<<3)
#define MD_RERUN     (1<<4)
#define MD_FFSCAN    (1<<6)
#define MD_SEPDVDL   (1<<7)

/* ROUTINES from md.c */
extern time_t do_md(FILE *log,t_commrec *cr,t_commrec *mcr,
		    int nfile,t_filenm fnm[],
		    bool bVerbose,bool bCompact,bool bVsites,
		    t_comm_vsites *vsitecomm,int stepout,
		    t_parm *parm,t_groups *grps,
		    t_topology *top,real ener[],t_fcdata *fcd,
		    t_state *state,rvec vold[],rvec vt[],rvec f[],
		    rvec buf[],t_mdatoms *mdatoms,
		    t_nsborder *nsb,t_nrnb nrnb[],
		    t_graph *graph,t_edsamyn *edyn,
		    t_forcerec *fr,rvec box_size,
		    int repl_ex_nst,int repl_ex_seed,
		    unsigned long Flags);

/* ROUTINES from minimize.c */
extern time_t do_steep(FILE *log,int nfile,t_filenm fnm[],
		       t_parm *parm,t_topology *top,
		       t_groups *grps,t_nsborder *nsb,
		       t_state *state,rvec grad[],rvec buf[],t_mdatoms *mdatoms,
		       tensor ekin,real ener[],t_fcdata *fcd,t_nrnb nrnb[],
		       bool bVerbose,bool bVsites,t_comm_vsites *vsitecomm,
		       t_commrec *cr,t_commrec *mcr,
		       t_graph *graph,t_forcerec *fr,rvec box_size);
/* Do steepest descents EM or something like that! */

extern time_t do_cg(FILE *log,int nfile,t_filenm fnm[],
		    t_parm *parm,t_topology *top,
		    t_groups *grps,t_nsborder *nsb,
		    t_state *state,rvec grad[],rvec buf[],t_mdatoms *mdatoms,
		    tensor ekin,real ener[],t_fcdata *fcd,t_nrnb nrnb[],
		    bool bVerbose,bool bVsites,t_comm_vsites *vsitecomm,
		    t_commrec *cr,t_commrec *mcr,
		    t_graph *graph,t_forcerec *fr,rvec box_size);
/* Do conjugate gradients EM! */

extern time_t do_lbfgs(FILE *log,int nfile,t_filenm fnm[],
		       t_parm *parm,t_topology *top,
		       t_groups *grps,t_nsborder *nsb, t_state *state,
		       rvec grad[],rvec buf[],t_mdatoms *mdatoms,
		       tensor ekin,real ener[],t_fcdata *fcd,t_nrnb nrnb[],
		       bool bVerbose,bool bVsites,t_comm_vsites *vsitecomm,
		       t_commrec *cr,t_commrec *mcr,
		       t_graph *graph,t_forcerec *fr,rvec box_size);
/* Do conjugate gradients EM! */


extern time_t do_nm(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
		    bool bVerbose,bool bCompact,int stepout,
		    t_parm *parm,t_groups *grps,
		    t_topology *top,real ener[],t_fcdata *fcd,
		    t_state *state,rvec vold[],rvec vt[],rvec f[],
		    rvec buf[],t_mdatoms *mdatoms,
		    t_nsborder *nsb,t_nrnb nrnb[],
		    t_graph *graph,t_edsamyn *edyn,
		    t_forcerec *fr,rvec box_size);
/* Do normal mode analysis */

/* ROUTINES from runner.c */
extern bool optRerunMDset (int nfile, t_filenm fnm[]);

extern void do_pbc_first(FILE *log,matrix box,rvec box_size,t_forcerec *fr,
			 t_graph *graph,rvec x[]);
		     
/* ROUTINES from stat.c */		
extern void global_stat(FILE *log,
			t_commrec *cr,real ener[],
			tensor fvir,tensor svir,
			t_grpopts *opts,t_groups *grps,
			t_nrnb *mynrnb,t_nrnb nrnb[],
			t_vcm *vcm,real *terminate);
/* Communicate statistics around the ring */

extern int write_traj(FILE *log,t_commrec *cr,char *traj,t_nsborder *nsb,
		      int step,real t,real lambda,t_nrnb nr_nb[],
		      int natoms,rvec *xx,rvec *vv,rvec *ff,matrix box);
/* Routine to output statusfiles during a run, as specified in
 * in parm->ir. If any of the pointers xx,vv,ff or ener is not NULL
 * it is written to the trajectory file.
 * Also write the energies etc. to the log file.
 * Returns the file handle (to be closed with close_trn).
 */

extern int do_per_step(int step,int nstep);
/* Return TRUE if io should be done */

extern int do_any_io(int step, t_inputrec *ir);

extern void write_xtc_traj(FILE *log,t_commrec *cr,
			   char *xtc_traj,t_nsborder *nsb,t_mdatoms *md,
			   int step,real t,rvec *xx,
			   matrix box,real prec);

extern void close_xtc_traj(void);

/* ROUTINES from sim_util.c */
extern void update_mdatoms(t_mdatoms *md,real lambda, bool bFirst);
/* Compute fields from mdatoms struct (invmass etc.) which may change
 * due to lambda dependent FEP calculations.
 * If bFirst all values are set, this is necessary once in the
 * first step.
 * You only have to call this routine again if lambda changes.
 */
 
extern void print_time(FILE *out,time_t start,int step,t_inputrec *ir);

extern time_t print_date_and_time(FILE *log,int pid,char *title);

extern void do_force(FILE *log,t_commrec *cr,t_commrec *mcr,
		     t_parm *parm,t_nsborder *nsb,
		     tensor vir_part,
		     int step,t_nrnb *nrnb,t_topology *top,t_groups *grps,
		     matrix box,rvec x[],rvec f[],rvec buf[],
		     t_mdatoms *mdatoms,real ener[],t_fcdata *fcd,
		     bool bVerbose,real lambda,t_graph *graph,
		     bool bNS,bool bNBFonly,t_forcerec *fr, rvec mu_tot,
		     bool bGatherOnly,real t,FILE *field);
extern void sum_lrforces(rvec f[],t_forcerec *fr,int start,int homenr);
		     
extern void calc_virial(FILE *log,int start,int homenr,rvec x[],rvec f[],
			tensor vir_part,tensor vir_el_recip,
			t_graph *graph,matrix box,
			t_nrnb *nrnb,const t_forcerec *fr);
			
extern void nstop_cm(FILE *log,t_commrec *cr,
		     int start,int nr_atoms,real mass[],rvec x[],rvec v[]);

extern void finish_run(FILE *log,t_commrec *cr,char *confout, t_nsborder *nsb,
		       t_topology *top, t_parm *parm,t_nrnb nrnb[],
		       double nodetime,double realtime,int step,
		       bool bWriteStat);

extern void calc_dispcorr(FILE *log,int eDispCorr,t_forcerec *fr,int natoms,
			  matrix box,tensor pres,tensor virial,real ener[]);
     

/* STUFF from init.c */
extern void write_parm(FILE *log,char *title,int pid,t_parm *parm);
/* Write parm for debugging */

typedef enum
{
  LIST_SCALARS	=0001,
  LIST_PARM	=0002,
  LIST_TOP	=0004,
  LIST_X	=0010,
  LIST_V	=0020,
  LIST_F	=0040,
  LIST_LOAD	=0100
} t_listitem;

extern void init_single(FILE *log,
                        t_parm *parm, char *tpbfile, t_topology *top,
			t_state *state,t_mdatoms **mdatoms,
			t_nsborder *nsb);
     /*
      * Allocates space for the topology (top), the coordinates x, the
      * velocities v, masses mass. Reads the parameters, topology,
      * coordinates and velocities from the file specified in tpbfile
      */

extern void distribute_parts(int left,int right,int pid,int nprocs,
                             t_parm *parm,char *tpbfile,int nstDlb);
     /*
      * Reads the parameters, topology, coordinates and velocities for the
      * multi processor version of the program from the file specified in
      * parm->files[STATUS_NM]. This file should also contain a so called
      * split descriptor which describes how to distribute particles over
      * the system. It then selects for all subsystems the appropriate data
      * and sends this to the processor using the left and right channels.
      * At last it sends its own subsystem down the ring where it is buffered.
      * Its own buffers for reading the data from the file are freed, and it
      * is now possible to reload this processor from the ring by using the
      * init_parts() routine.
      * The routine also creates a renum array which can be used for writing
      * out the x,v and f for analysis purpose.
      */

extern void init_parts(FILE *log,t_commrec *cr,
		       t_parm *parm,t_topology *top,
		       t_state *state,t_mdatoms **mdatoms,
		       t_nsborder *nsb,int list,
		       bool *bParallelVsites,
		       t_comm_vsites *vsitecomm);
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

extern void do_shakefirst(FILE *log,bool bTYZ,real ener[],
			  t_parm *parm,t_nsborder *nsb,t_mdatoms *md,
			  t_state *state,rvec vold[],rvec buf[],rvec f[],
			  t_graph *graph,t_commrec *cr,t_nrnb *nrnb,
			  t_groups *grps,t_forcerec *fr,t_topology *top,
			  t_edsamyn *edyn,t_pull *pulldata);
			  
extern void dynamic_load_balancing(bool bVerbose,t_commrec *cr,real capacity[],
				   int dimension,t_mdatoms *md,t_topology *top,
				   rvec x[],rvec v[],matrix box);
/* Perform load balancing, i.e. split the particles over processors
 * based on their coordinates in the "dimension" direction.
 */
				   
extern void mdrunner(t_commrec *cr,t_commrec *mcr,int nfile,t_filenm fnm[],
		     bool bVerbose,bool bCompact,
		     int nDlb,int nstepout,
		     t_edsamyn *edyn,int repl_ex_nst,int repl_ex_seed,
		     unsigned long Flags);
/* Driver routine, that calls the different methods */

extern void init_md(t_commrec *cr,t_inputrec *ir,tensor box,real *t,real *t0,
		    real *lambda,real *lam0,
		    t_nrnb *mynrnb,bool *bTYZ,t_topology *top,
		    int nfile,t_filenm fnm[],char **traj,
		    char **xtc_traj,int *fp_ene,
		    FILE **fp_dgdl,FILE **fp_field,
		    t_mdebin **mdebin,t_groups *grps,
		    tensor force_vir,tensor shake_vir,
		    t_mdatoms *mdatoms,rvec mu_tot,
		    bool *bNEMD,bool *bSimAnn,t_vcm **vcm,t_nsborder *nsb);
/* Routine in sim_util.c */

extern void init_em(FILE *log,const char *title,t_parm *parm,
		    real *lambda,t_nrnb *mynrnb,rvec mu_tot,
		    matrix box,rvec box_size,
		    t_forcerec *fr,t_mdatoms *mdatoms,t_topology *top,
		    t_nsborder *nsb,
		    t_commrec *cr,t_vcm **vcm,int *start,int *end);
/* Routine in minimize.c */
     
#endif	/* _mdrun_h */
