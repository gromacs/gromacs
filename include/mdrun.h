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

#ifndef _mdrun_h
#define _mdrun_h

static char *SRCID_mdrun_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) do_md.h 1.12 03 Mar 1996"
#endif /* HAVE_IDENT */
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

/* ROUTINES from md.c */
extern time_t do_md(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
		    bool bVerbose,bool bCompact,bool bDummies,int stepout,
		    t_parm *parm,t_groups *grps,
		    t_topology *top,real ener[],
		    rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
		    rvec buf[],t_mdatoms *mdatoms,
		    t_nsborder *nsb,t_nrnb nrnb[],
		    t_graph *graph,t_edsamyn *edyn,
		    t_forcerec *fr,rvec box_size);

/* ROUTINES from nm.c */
extern time_t do_nm(FILE *log,t_commrec *cr,int nfile,t_filenm fnm[],
		    bool bVerbose,bool bCompact,int stepout,
		    t_parm *parm,t_groups *grps,
		    t_topology *top,real ener[],
		    rvec x[],rvec vold[],rvec v[],rvec vt[],rvec f[],
		    rvec buf[],t_mdatoms *mdatoms,
		    t_nsborder *nsb,t_nrnb nrnb[],
		    t_graph *graph,t_edsamyn *edyn,
		    t_forcerec *fr,rvec box_size);

/* ROUTINES from steep.c */
extern real f_norm(int left,int right,int nprocs,
		   int start,int end,rvec grad[]);
/* Calculates norm of force */

extern real f_max(int left,int right,int nprocs,
		  int start,int end,rvec grad[]);
/* Calculates max force */

extern time_t do_steep(FILE *log,int nfile,t_filenm fnm[],
		       t_parm *parm,t_topology *top,
		       t_groups *grps,t_nsborder *nsb,
		       rvec x[],rvec grad[],rvec buf[],t_mdatoms *mdatoms,
		       tensor ekin,real ener[],t_nrnb nrnb[],
		       bool bVerbose,bool bDummies,t_commrec *cr,
		       t_graph *graph,t_forcerec *fr,rvec box_size);
/* Do steepest descents EM or something like that! */

/* ROUTINES from congrad.c */
extern time_t do_cg(FILE *log,int nfile,t_filenm fnm[],
		    t_parm *parm,t_topology *top,
		    t_groups *grps,t_nsborder *nsb,
		    rvec x[],rvec grad[],rvec buf[],t_mdatoms *mdatoms,
		    tensor ekin,real ener[],t_nrnb nrnb[],
		    bool bVerbose,bool bDummies,
		    t_commrec *cr,t_graph *graph,t_forcerec *fr,
		    rvec box_size);
/* Do conjugate gradients EM! */

/* ROUTINES from runner.c */
extern bool optRerunMDset (int nfile, t_filenm fnm[]);

extern void mdrunner(t_commrec *cr,int nfile,t_filenm fnm[],bool bVerbose,
		     bool bCompact,int nDlb,bool bNM,int nstepout,
		     t_edsamyn *edyn);
		    
/* Initialization routines to make maintainance easier */ 
extern void init_md(t_commrec *cr,
		    t_inputrec *ir,real *t,real *t0,real *lambda,real *lam0,
		    real *SAfactor,t_nrnb *mynrnb,bool *bTYZ,t_topology *top,
		    int nfile,t_filenm fnm[],char **traj,char **xtc_traj,
		    int *fp_ene,FILE **fp_dgdl,
		    t_mdebin **mdebin,t_groups *grps,rvec vcm,tensor force_vir,
		    tensor shake_vir,t_mdatoms *mdatoms);
		    
extern void do_pbc_first(FILE *log,t_parm *parm,rvec box_size,t_forcerec *fr,
			 t_graph *graph,rvec x[]);
		     
extern void get_cmparm(t_inputrec *ir,int step,bool *bStopCM,bool *bStopRot);
/* Initiate center of mass removal parameters */

void set_pot_bools(t_inputrec *ir,t_topology *top,
		   bool *bLR,bool *bLJLR,bool *bBHAM,bool *b14);
/* Initiate some bools for the potential energy calculation */

/* ROUTINES from stat.c */		
extern void global_stat(FILE *log,
			t_commrec *cr,real ener[],
			tensor fvir,tensor svir,
			t_grpopts *opts,t_groups *grps,
			t_nrnb *mynrnb,t_nrnb nrnb[],
			rvec vcm,rvec mu_tot,real *terminate);
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
extern void init_mdatoms(t_mdatoms *md,real lambda,bool bFirst);
/* Compute fields from mdatoms struct (invmass etc.) which may change
 * due to lambda dependent FEP calculations.
 * If bFirst all values are set, this is necessary once in the
 * first step.
 */
 
extern void print_time(FILE *out,time_t start,int step,t_inputrec *ir);

extern time_t print_date_and_time(FILE *log,int pid,char *title);

extern void do_force(FILE *log,t_commrec *cr,
		     t_parm *parm,t_nsborder *nsb,tensor vir_part,
		     int step,t_nrnb *nrnb,t_topology *top,t_groups *grps,
		     rvec x[],rvec v[],rvec f[],rvec buf[],
		     t_mdatoms *mdatoms,real ener[],bool bVerbose,
		     real lambda,t_graph *graph,
		     bool bNS,bool bNBFonly,t_forcerec *fr);

extern void nstop_cm(FILE *log,t_commrec *cr,
		     int start,int nr_atoms,real mass[],rvec x[],rvec v[]);

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
			rvec **x,rvec **v,t_mdatoms **mdatoms,
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
		       rvec **x,rvec **v,t_mdatoms **mdatoms,
		       t_nsborder *nsb,int list);
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
 
extern double cpu_time(void);
/* Return the cpu time so far in seconds. */

extern void do_shakefirst(FILE *log,bool bTYZ,real lambda,real ener[],
			  t_parm *parm,t_nsborder *nsb,t_mdatoms *md,
			  rvec x[],rvec vold[],rvec buf[],rvec f[],
			  rvec v[],t_graph *graph,t_commrec *cr,t_nrnb *nrnb,
			  t_groups *grps,t_forcerec *fr,t_topology *top,
			  t_edsamyn *edyn,t_pull *pulldata);
			  
extern void get_cmparm(t_inputrec *ir,int step,bool *bStopCM,bool *bStopRot);
/* Determine from the input whether or not to stop center of mass motion */

#endif	/* _mdrun_h */
