/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef _mdrun_h
#define _mdrun_h

#include <stdio.h>
#include <time.h>
#include "visibility.h"
#include "typedefs.h"
#include "network.h"
#include "sim_util.h"
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
#include "types/membedt.h"
#include "types/globsig.h"


#ifdef GMX_THREAD_MPI
#include "thread_mpi/threads.h"
#endif

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
#define MD_APPENDFILESSET (1<<21)
#define MD_KEEPANDNUMCPT  (1<<16)
#define MD_READ_EKIN      (1<<17)
#define MD_STARTFROMCPT   (1<<18)
#define MD_RESETCOUNTERSHALFWAY (1<<19)
#define MD_TUNEPME        (1<<20)
#define MD_TESTVERLET     (1<<22)

/* The options for the domain decomposition MPI task ordering */
enum {
    ddnoSEL, ddnoINTERLEAVE, ddnoPP_PME, ddnoCARTESIAN, ddnoNR
};

/* Variables for temporary use with the deform option,
 * used in runner.c and md.c.
 * (These variables should be stored in the tpx file.)
 */
extern gmx_large_int_t     deform_init_init_step_tpx;
extern matrix              deform_init_box_tpx;
#ifdef GMX_THREAD_MPI
extern tMPI_Thread_mutex_t deform_init_box_mutex;

/* The minimum number of atoms per tMPI thread. With fewer atoms than this,
 * the number of threads will get lowered.
 */
#define MIN_ATOMS_PER_MPI_THREAD    90
#define MIN_ATOMS_PER_GPU           900
#endif


typedef double gmx_integrator_t (FILE *log, t_commrec *cr,
                                 int nfile, const t_filenm fnm[],
                                 const output_env_t oenv, gmx_bool bVerbose,
                                 gmx_bool bCompact, int nstglobalcomm,
                                 gmx_vsite_t *vsite, gmx_constr_t constr,
                                 int stepout,
                                 t_inputrec *inputrec,
                                 gmx_mtop_t *mtop, t_fcdata *fcd,
                                 t_state *state,
                                 t_mdatoms *mdatoms,
                                 t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                                 gmx_edsam_t ed,
                                 t_forcerec *fr,
                                 int repl_ex_nst, int repl_ex_nex, int repl_ex_seed,
                                 gmx_membed_t membed,
                                 real cpt_period, real max_hours,
                                 const char *deviceOptions,
                                 unsigned long Flags,
                                 gmx_runtime_t *runtime);

/* ROUTINES from md.c */

gmx_integrator_t do_md;


/* ROUTINES from minimize.c */

GMX_LIBMD_EXPORT
gmx_integrator_t do_steep;
/* Do steepest descents EM */

GMX_LIBMD_EXPORT
gmx_integrator_t do_cg;
/* Do conjugate gradient EM */

GMX_LIBMD_EXPORT
gmx_integrator_t do_lbfgs;
/* Do conjugate gradient L-BFGS */

GMX_LIBMD_EXPORT
gmx_integrator_t do_nm;
/* Do normal mode analysis */

/* ROUTINES from tpi.c */

GMX_LIBMD_EXPORT
gmx_integrator_t do_tpi;
/* Do test particle insertion */

void init_npt_masses(t_inputrec *ir, t_state *state, t_extmass *MassQ, gmx_bool bInit);

GMX_LIBMD_EXPORT
void init_expanded_ensemble(gmx_bool bStateFromCP, t_inputrec *ir, gmx_rng_t *mcrng, df_history_t *dfhist);

GMX_LIBMD_EXPORT
int ExpandedEnsembleDynamics(FILE *log, t_inputrec *ir, gmx_enerdata_t *enerd,
                             t_state *state, t_extmass *MassQ, int fep_state, df_history_t *dfhist,
                             gmx_large_int_t step, gmx_rng_t mcrng,
                             rvec *v, t_mdatoms *mdatoms);

GMX_LIBMD_EXPORT
void PrintFreeEnergyInfoToFile(FILE *outfile, t_lambda *fep, t_expanded *expand, t_simtemp *simtemp, df_history_t *dfhist,
                               int fep_state, int frequency, gmx_large_int_t step);

GMX_LIBMD_EXPORT
void get_mc_state(gmx_rng_t rng, t_state *state);

GMX_LIBMD_EXPORT
void set_mc_state(gmx_rng_t rng, t_state *state);

/* check the version */
GMX_LIBMD_EXPORT
void check_ir_old_tpx_versions(t_commrec *cr, FILE *fplog,
                               t_inputrec *ir, gmx_mtop_t *mtop);

/* Allocate and initialize node-local state entries. */
GMX_LIBMD_EXPORT
void set_state_entries(t_state *state, const t_inputrec *ir, int nnodes);

/* Broadcast the data for a simulation, and allocate node-specific settings
   such as rng generators. */
GMX_LIBMD_EXPORT
void init_parallel(FILE *log, t_commrec *cr, t_inputrec *inputrec,
                   gmx_mtop_t *mtop);

int mdrunner(gmx_hw_opt_t *hw_opt,
             FILE *fplog, t_commrec *cr, int nfile,
             const t_filenm fnm[], const output_env_t oenv, gmx_bool bVerbose,
             gmx_bool bCompact, int nstglobalcomm, ivec ddxyz, int dd_node_order,
             real rdd, real rconstr, const char *dddlb_opt, real dlb_scale,
             const char *ddcsx, const char *ddcsy, const char *ddcsz,
             const char *nbpu_opt,
             gmx_large_int_t nsteps_cmdline, int nstepout, int resetstep,
             int nmultisim, int repl_ex_nst, int repl_ex_nex,
             int repl_ex_seed, real pforce, real cpt_period, real max_hours,
             const char *deviceOptions, unsigned long Flags);
/* Driver routine, that calls the different methods */

#ifdef __cplusplus
}
#endif

#endif  /* _mdrun_h */
