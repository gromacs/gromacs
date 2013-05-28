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

#ifndef _sim_util_h
#define _sim_util_h

#include <time.h>
#include "visibility.h"
#include "typedefs.h"
#include "enxio.h"
#include "mdebin.h"
#include "update.h"
#include "vcm.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    t_fileio   *fp_trn;
    t_fileio   *fp_xtc;
    int         xtc_prec;
    ener_file_t fp_ene;
    const char *fn_cpt;
    gmx_bool    bKeepAndNumCPT;
    int         eIntegrator;
    gmx_bool    bExpanded;
    int         elamstats;
    int         simulation_part;
    FILE       *fp_dhdl;
    FILE       *fp_field;
} gmx_mdoutf_t;

typedef struct gmx_global_stat *gmx_global_stat_t;

typedef struct {
    double          real;
#ifdef GMX_CRAY_XT3
    double          proc;
#else
    clock_t         proc;
#endif
    double          realtime;
    double          proctime;
    double          time_per_step;
    double          last;
    gmx_large_int_t nsteps_done;
} gmx_runtime_t;


GMX_LIBMD_EXPORT
void do_pbc_first(FILE *log, matrix box, t_forcerec *fr,
                  t_graph *graph, rvec x[]);

GMX_LIBMD_EXPORT
void do_pbc_first_mtop(FILE *fplog, int ePBC, matrix box,
                       gmx_mtop_t *mtop, rvec x[]);

GMX_LIBMD_EXPORT
void do_pbc_mtop(FILE *fplog, int ePBC, matrix box,
                 gmx_mtop_t *mtop, rvec x[]);



/* ROUTINES from stat.c */
GMX_LIBMD_EXPORT
gmx_global_stat_t global_stat_init(t_inputrec *ir);

void global_stat_destroy(gmx_global_stat_t gs);

void global_stat(FILE *log, gmx_global_stat_t gs,
                 t_commrec *cr, gmx_enerdata_t *enerd,
                 tensor fvir, tensor svir, rvec mu_tot,
                 t_inputrec *inputrec,
                 gmx_ekindata_t *ekind,
                 gmx_constr_t constr, t_vcm *vcm,
                 int nsig, real *sig,
                 gmx_mtop_t *top_global, t_state *state_local,
                 gmx_bool bSumEkinhOld, int flags);
/* Communicate statistics over cr->mpi_comm_mysim */

gmx_mdoutf_t *init_mdoutf(int nfile, const t_filenm fnm[],
                          int mdrun_flags,
                          const t_commrec *cr, const t_inputrec *ir,
                          const output_env_t oenv);
/* Returns a pointer to a data structure with all output file pointers
 * and names required by mdrun.
 */

GMX_LIBMD_EXPORT
void done_mdoutf(gmx_mdoutf_t *of);
/* Close all open output files and free the of pointer */

#define MDOF_X   (1<<0)
#define MDOF_V   (1<<1)
#define MDOF_F   (1<<2)
#define MDOF_XTC (1<<3)
#define MDOF_CPT (1<<4)

GMX_LIBMD_EXPORT
void write_traj(FILE *fplog, t_commrec *cr,
                gmx_mdoutf_t *of,
                int mdof_flags,
                gmx_mtop_t *top_global,
                gmx_large_int_t step, double t,
                t_state *state_local, t_state *state_global,
                rvec *f_local, rvec *f_global,
                int *n_xtc, rvec **x_xtc);
/* Routine that writes frames to trn, xtc and/or checkpoint.
 * What is written is determined by the mdof_flags defined above.
 * Data is collected to the master node only when necessary.
 */

GMX_LIBMD_EXPORT
int do_per_step(gmx_large_int_t step, gmx_large_int_t nstep);
/* Return TRUE if io should be done */

/* ROUTINES from sim_util.c */

GMX_LIBMD_EXPORT
double gmx_gettime();

GMX_LIBMD_EXPORT
void print_time(FILE *out, gmx_runtime_t *runtime,
                gmx_large_int_t step, t_inputrec *ir, t_commrec *cr);

GMX_LIBMD_EXPORT
void runtime_start(gmx_runtime_t *runtime);

GMX_LIBMD_EXPORT
void runtime_end(gmx_runtime_t *runtime);

GMX_LIBMD_EXPORT
void runtime_upd_proc(gmx_runtime_t *runtime);
/* The processor time should be updated every once in a while,
 * since on 32-bit manchines it loops after 72 minutes.
 */

GMX_LIBMD_EXPORT
void print_date_and_time(FILE *log, int pid, const char *title,
                         const gmx_runtime_t *runtime);

GMX_LIBMD_EXPORT
void print_start(FILE *fplog, t_commrec *cr, gmx_runtime_t *runtime,
                 const char *name);

GMX_LIBMD_EXPORT
void finish_run(FILE *log, t_commrec *cr, const char *confout,
                t_inputrec *inputrec,
                t_nrnb nrnb[], gmx_wallcycle_t wcycle,
                gmx_runtime_t *runtime,
                wallclock_gpu_t *gputimes,
                int omp_nth_pp,
                gmx_bool bWriteStat);

void calc_enervirdiff(FILE *fplog, int eDispCorr, t_forcerec *fr);

void calc_dispcorr(FILE *fplog, t_inputrec *ir, t_forcerec *fr,
                   gmx_large_int_t step, int natoms,
                   matrix box, real lambda, tensor pres, tensor virial,
                   real *prescorr, real *enercorr, real *dvdlcorr);

void initialize_lambdas(FILE *fplog, t_inputrec *ir, int *fep_state, real *lambda, double *lam0);

GMX_LIBMD_EXPORT
void do_constrain_first(FILE *log, gmx_constr_t constr,
                        t_inputrec *inputrec, t_mdatoms *md,
                        t_state *state, rvec *f,
                        t_graph *graph, t_commrec *cr, t_nrnb *nrnb,
                        t_forcerec *fr, gmx_localtop_t *top, tensor shake_vir);

GMX_LIBMD_EXPORT
void init_md(FILE *fplog,
             t_commrec *cr, t_inputrec *ir, const output_env_t oenv,
             double *t, double *t0,
             real *lambda, int *fep_state, double *lam0,
             t_nrnb *nrnb, gmx_mtop_t *mtop,
             gmx_update_t *upd,
             int nfile, const t_filenm fnm[],
             gmx_mdoutf_t **outf, t_mdebin **mdebin,
             tensor force_vir, tensor shake_vir,
             rvec mu_tot,
             gmx_bool *bSimAnn, t_vcm **vcm,
             t_state *state, unsigned long Flags);
/* Routine in sim_util.c */

#ifdef __cplusplus
}
#endif

#endif  /* _sim_util_h */
