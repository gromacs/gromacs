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

#ifndef _mdrun_h
#define _mdrun_h

#include <stdio.h>
#include <time.h>
#include "typedefs.h"
#include "network.h"
#include "tgroup.h"
#include "filenm.h"
#include "mshift.h"
#include "force.h"
#include "edsam.h"
#include "mdebin.h"
#include "vcm.h"
#include "vsite.h"
#include "pull.h"
#include "update.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MD_POLARISE       (1<<2)
#define MD_IONIZE         (1<<3)
#define MD_RERUN          (1<<4)
#define MD_RERUN_VSITE    (1<<5)
#define MD_FFSCAN         (1<<6)
#define MD_SEPPOT         (1<<7)
#define MD_PARTDEC        (1<<9)
#define MD_DDBONDCHECK    (1<<10)
#define MD_DDBONDCOMM     (1<<11)
#define MD_CONFOUT        (1<<12)
#define MD_REPRODUCIBLE   (1<<13)
#define MD_READ_RNG       (1<<14)
#define MD_APPENDFILES    (1<<15)
#define MD_KEEPANDNUMCPT  (1<<16)
#define MD_READ_EKIN      (1<<17)
#define MD_STARTFROMCPT   (1<<18)
#define MD_RESETCOUNTERSHALFWAY (1<<19)

/* Define a number of flags to better control the information
 * passed to compute_globals in md.c and global_stat.
 */

/* We are rerunning the simulation */
#define CGLO_RERUNMD        (1<<1)
/* we are computing the kinetic energy from average velocities */
#define CGLO_EKINAVEVEL     (1<<2)
/* we are removing the center of mass momenta */
#define CGLO_STOPCM         (1<<3)
/* bGStat is defined in do_md */
#define CGLO_GSTAT          (1<<4)
/* Sum the energy terms in global computation */
#define CGLO_ENERGY         (1<<6)
/* Sum the kinetic energy terms in global computation */
#define CGLO_TEMPERATURE    (1<<7)
/* Sum the kinetic energy terms in global computation */
#define CGLO_PRESSURE       (1<<8)
/* Sum the constraint term in global computation */
#define CGLO_CONSTRAINT     (1<<9)
/* we are using an integrator that requires iteration over some steps - currently not used*/
#define CGLO_ITERATE        (1<<10)
/* it is the first time we are iterating (or, only once through is required */
#define CGLO_FIRSTITERATE   (1<<11)
/* Reading ekin from the trajectory */
#define CGLO_READEKIN       (1<<12)
/* we need to reset the ekin rescaling factor here */
#define CGLO_SCALEEKIN      (1<<13)
  
enum {
  ddnoSEL, ddnoINTERLEAVE, ddnoPP_PME, ddnoCARTESIAN, ddnoNR
};

typedef struct {
  double real;
#ifdef GMX_CRAY_XT3
  double proc;
#else
  clock_t proc;
#endif
  double realtime;
  double proctime;
  double time_per_step;
  double last;
  gmx_large_int_t nsteps_done;
} gmx_runtime_t;

typedef struct {
  t_fileio *fp_trn;
  t_fileio *fp_xtc;
  int  xtc_prec;
  ener_file_t fp_ene;
  const char *fn_cpt;
  gmx_bool bKeepAndNumCPT;
  int  eIntegrator;
  int  simulation_part;
  FILE *fp_dhdl;
  FILE *fp_field;
} gmx_mdoutf_t;

/* Variables for temporary use with the deform option,
 * used in runner.c and md.c.
 * (These variables should be stored in the tpx file.)
 */
extern gmx_large_int_t     deform_init_init_step_tpx;
extern matrix              deform_init_box_tpx;
#ifdef GMX_THREADS
extern tMPI_Thread_mutex_t deform_init_box_mutex;

/* The minimum number of atoms per thread. With fewer atoms than this,
 * the number of threads will get lowered.
 */
#define MIN_ATOMS_PER_THREAD    90
#endif


typedef double gmx_integrator_t(FILE *log,t_commrec *cr,
				int nfile,const t_filenm fnm[],
				const output_env_t oenv, gmx_bool bVerbose,
                                gmx_bool bCompact, int nstglobalcomm,
				gmx_vsite_t *vsite,gmx_constr_t constr,
				int stepout,
				t_inputrec *inputrec,
				gmx_mtop_t *mtop,t_fcdata *fcd,
				t_state *state,
				t_mdatoms *mdatoms,
				t_nrnb *nrnb,gmx_wallcycle_t wcycle,
				gmx_edsam_t ed, 
				t_forcerec *fr,
				int repl_ex_nst,int repl_ex_seed,
				real cpt_period,real max_hours,
				const char *deviceOptions,
				unsigned long Flags,
				gmx_runtime_t *runtime);

typedef struct gmx_global_stat *gmx_global_stat_t;

/* ROUTINES from md.c */

gmx_integrator_t do_md;

gmx_integrator_t do_md_openmm;



/* ROUTINES from minimize.c */

gmx_integrator_t do_steep;
/* Do steepest descents EM */

gmx_integrator_t do_cg;
/* Do conjugate gradient EM */

gmx_integrator_t do_lbfgs;
/* Do conjugate gradient L-BFGS */

gmx_integrator_t do_nm;
/* Do normal mode analysis */

/* ROUTINES from tpi.c */

gmx_integrator_t do_tpi;
/* Do test particle insertion */


/* ROUTINES from md_support.c */

/* return the number of steps between global communcations */
int check_nstglobalcomm(FILE *fplog,t_commrec *cr,
                        int nstglobalcomm,t_inputrec *ir);

/* check whether an 'nst'-style parameter p is a multiple of nst, and
   set it to be one if not, with a warning. */
void check_nst_param(FILE *fplog,t_commrec *cr,
                     const char *desc_nst,int nst,
                     const char *desc_p,int *p);

/* check which of the multisim simulations has the shortest number of
   steps and return that number of nsteps */
gmx_large_int_t get_multisim_nsteps(const t_commrec *cr,
                                    gmx_large_int_t nsteps);

void rerun_parallel_comm(t_commrec *cr,t_trxframe *fr,
                         gmx_bool *bNotLastFrame);

/* get the conserved energy associated with the ensemble type*/
real compute_conserved_from_auxiliary(t_inputrec *ir, t_state *state,           
                                      t_extmass *MassQ);

/* reset all cycle and time counters. */
void reset_all_counters(FILE *fplog,t_commrec *cr,
                        gmx_large_int_t step,
                        gmx_large_int_t *step_rel,t_inputrec *ir,
                        gmx_wallcycle_t wcycle,t_nrnb *nrnb,
                        gmx_runtime_t *runtime);



/* ROUTINES from sim_util.c */
void do_pbc_first(FILE *log,matrix box,t_forcerec *fr,
			 t_graph *graph,rvec x[]);

void do_pbc_first_mtop(FILE *fplog,int ePBC,matrix box,
			      gmx_mtop_t *mtop,rvec x[]);

void do_pbc_mtop(FILE *fplog,int ePBC,matrix box,
			gmx_mtop_t *mtop,rvec x[]);


		     
/* ROUTINES from stat.c */
gmx_global_stat_t global_stat_init(t_inputrec *ir);

void global_stat_destroy(gmx_global_stat_t gs);

void global_stat(FILE *log,gmx_global_stat_t gs,
			t_commrec *cr,gmx_enerdata_t *enerd,
			tensor fvir,tensor svir,rvec mu_tot,
			t_inputrec *inputrec,
			gmx_ekindata_t *ekind,
			gmx_constr_t constr,t_vcm *vcm,
			int nsig,real *sig,
			gmx_mtop_t *top_global, t_state *state_local, 
			gmx_bool bSumEkinhOld, int flags);
/* Communicate statistics over cr->mpi_comm_mysim */

gmx_mdoutf_t *init_mdoutf(int nfile,const t_filenm fnm[],
				 int mdrun_flags,
				 const t_commrec *cr,const t_inputrec *ir,
				 const output_env_t oenv);
/* Returns a pointer to a data structure with all output file pointers
 * and names required by mdrun.
 */

void done_mdoutf(gmx_mdoutf_t *of);
/* Close all open output files and free the of pointer */

#define MDOF_X   (1<<0)
#define MDOF_V   (1<<1)
#define MDOF_F   (1<<2)
#define MDOF_XTC (1<<3)
#define MDOF_CPT (1<<4)

void write_traj(FILE *fplog,t_commrec *cr,
		       gmx_mdoutf_t *of,
		       int mdof_flags,
		       gmx_mtop_t *top_global,
		       gmx_large_int_t step,double t,
		       t_state *state_local,t_state *state_global,
		       rvec *f_local,rvec *f_global,
		       int *n_xtc,rvec **x_xtc);
/* Routine that writes frames to trn, xtc and/or checkpoint.
 * What is written is determined by the mdof_flags defined above.
 * Data is collected to the master node only when necessary.
 */

int do_per_step(gmx_large_int_t step,gmx_large_int_t nstep);
/* Return TRUE if io should be done */

int do_any_io(int step, t_inputrec *ir);

/* ROUTINES from sim_util.c */

double gmx_gettime();

void print_time(FILE *out, gmx_runtime_t *runtime,
                       gmx_large_int_t step,t_inputrec *ir, t_commrec *cr);

void runtime_start(gmx_runtime_t *runtime);

void runtime_end(gmx_runtime_t *runtime);

void runtime_upd_proc(gmx_runtime_t *runtime);
/* The processor time should be updated every once in a while,
 * since on 32-bit manchines it loops after 72 minutes.
 */
  
void print_date_and_time(FILE *log,int pid,const char *title,
				const gmx_runtime_t *runtime);
  
void nstop_cm(FILE *log,t_commrec *cr,
		     int start,int nr_atoms,real mass[],rvec x[],rvec v[]);

void finish_run(FILE *log,t_commrec *cr,const char *confout,
		       t_inputrec *inputrec,
		       t_nrnb nrnb[],gmx_wallcycle_t wcycle,
		       gmx_runtime_t *runtime,
		       gmx_bool bWriteStat);

void calc_enervirdiff(FILE *fplog,int eDispCorr,t_forcerec *fr);

void calc_dispcorr(FILE *fplog,t_inputrec *ir,t_forcerec *fr,
			  gmx_large_int_t step, int natoms, 
			  matrix box,real lambda,tensor pres,tensor virial,
			  real *prescorr, real *enercorr, real *dvdlcorr);

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

void check_nnodes_top(char *fn,t_topology *top);
/* Reset the tpr file to work with one node if necessary */


/* check the version */
void check_ir_old_tpx_versions(t_commrec *cr,FILE *fplog,
                               t_inputrec *ir,gmx_mtop_t *mtop);

/* Allocate and initialize node-local state entries. */
void set_state_entries(t_state *state,const t_inputrec *ir,int nnodes);

/* Broadcast the data for a simulation, and allocate node-specific settings
   such as rng generators. */
void init_parallel(FILE *log, t_commrec *cr, t_inputrec *inputrec,
                          gmx_mtop_t *mtop);


void do_constrain_first(FILE *log,gmx_constr_t constr,
			       t_inputrec *inputrec,t_mdatoms *md,
			       t_state *state,rvec *f,
			       t_graph *graph,t_commrec *cr,t_nrnb *nrnb,
			       t_forcerec *fr, gmx_localtop_t *top, tensor shake_vir); 
			  
void dynamic_load_balancing(gmx_bool bVerbose,t_commrec *cr,real capacity[],
				   int dimension,t_mdatoms *md,t_topology *top,
				   rvec x[],rvec v[],matrix box);
/* Perform load balancing, i.e. split the particles over processors
 * based on their coordinates in the "dimension" direction.
 */

int multisim_min(const gmx_multisim_t *ms,int nmin,int n);
/* Set an appropriate value for n across the whole multi-simulation */

int multisim_nstsimsync(const t_commrec *cr,
			const t_inputrec *ir,int repl_ex_nst);
/* Determine the interval for inter-simulation communication */
				   
void init_global_signals(globsig_t *gs,const t_commrec *cr,
			 const t_inputrec *ir,int repl_ex_nst);
/* Constructor for globsig_t */

void copy_coupling_state(t_state *statea,t_state *stateb,
			 gmx_ekindata_t *ekinda,gmx_ekindata_t *ekindb, t_grpopts* opts);
/* Copy stuff from state A to state B */

void compute_globals(FILE *fplog, gmx_global_stat_t gstat, t_commrec *cr, t_inputrec *ir,
		     t_forcerec *fr, gmx_ekindata_t *ekind,
		     t_state *state, t_state *state_global, t_mdatoms *mdatoms,
		     t_nrnb *nrnb, t_vcm *vcm, gmx_wallcycle_t wcycle,
		     gmx_enerdata_t *enerd,tensor force_vir, tensor shake_vir, tensor total_vir,
		     tensor pres, rvec mu_tot, gmx_constr_t constr,
		     globsig_t *gs,gmx_bool bInterSimGS,
		     matrix box, gmx_mtop_t *top_global, real *pcurr,
		     int natoms, gmx_bool *bSumEkinhOld, int flags);
/* Compute global variables during integration */

int mdrunner(int nthreads_requested, FILE *fplog,t_commrec *cr,int nfile,
             const t_filenm fnm[], const output_env_t oenv, gmx_bool bVerbose,
             gmx_bool bCompact, int nstglobalcomm, ivec ddxyz,int dd_node_order,
             real rdd, real rconstr, const char *dddlb_opt,real dlb_scale,
	     const char *ddcsx,const char *ddcsy,const char *ddcsz,
	     int nstepout, int resetstep, int nmultisim, int repl_ex_nst,
             int repl_ex_seed, real pforce,real cpt_period,real max_hours,
	     const char *deviceOptions, unsigned long Flags);
/* Driver routine, that calls the different methods */

void md_print_warning(const t_commrec *cr,FILE *fplog,const char *buf);
/* Print a warning message to stderr on the master node
 * and to fplog if fplog!=NULL.
 */

void init_md(FILE *fplog,
		    t_commrec *cr,t_inputrec *ir, const output_env_t oenv, 
		    double *t,double *t0,
		    real *lambda,double *lam0,
		    t_nrnb *nrnb,gmx_mtop_t *mtop,
		    gmx_update_t *upd,
		    int nfile,const t_filenm fnm[],
		    gmx_mdoutf_t **outf,t_mdebin **mdebin,
		    tensor force_vir,tensor shake_vir,
		    rvec mu_tot,
		    gmx_bool *bSimAnn,t_vcm **vcm, 
		    t_state *state, unsigned long Flags);
  /* Routine in sim_util.c */

#ifdef __cplusplus
}
#endif

#endif	/* _mdrun_h */
