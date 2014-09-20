/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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

#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/mdoutf.h"
#include "gromacs/legacyheaders/mdebin.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/update.h"
#include "gromacs/legacyheaders/vcm.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"

#ifdef __cplusplus
extern "C" {
#endif

struct t_graph;

typedef struct gmx_global_stat *gmx_global_stat_t;

void do_pbc_first(FILE *log, matrix box, t_forcerec *fr,
                  struct t_graph *graph, rvec x[]);

void do_pbc_first_mtop(FILE *fplog, int ePBC, matrix box,
                       gmx_mtop_t *mtop, rvec x[]);

void do_pbc_mtop(FILE *fplog, int ePBC, matrix box,
                 gmx_mtop_t *mtop, rvec x[]);



/* ROUTINES from stat.c */
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

int do_per_step(gmx_int64_t step, gmx_int64_t nstep);
/* Return TRUE if io should be done */

/* ROUTINES from sim_util.c */

void print_time(FILE *out, gmx_walltime_accounting_t walltime_accounting,
                gmx_int64_t step, t_inputrec *ir, t_commrec *cr);

/*! \brief Print date, time, MPI rank and a description of this point
 * in time.
 *
 * \param[in] log       logfile, or NULL to suppress output
 * \param[in] rank      MPI rank to include in the output
 * \param[in] title     Description to include in the output
 * \param[in] the_time  Seconds since the epoch, e.g. as reported by gmx_gettime
 */
void print_date_and_time(FILE *log, int rank, const char *title,
                         double the_time);

void print_start(FILE *fplog, t_commrec *cr,
                 gmx_walltime_accounting_t walltime_accounting,
                 const char *name);

void finish_run(FILE *log, t_commrec *cr,
                t_inputrec *inputrec,
                t_nrnb nrnb[], gmx_wallcycle_t wcycle,
                gmx_walltime_accounting_t walltime_accounting,
                struct nonbonded_verlet_t *nbv,
                gmx_bool bWriteStat);

void calc_enervirdiff(FILE *fplog, int eDispCorr, t_forcerec *fr);

void calc_dispcorr(t_inputrec *ir, t_forcerec *fr,
                   int natoms,
                   matrix box, real lambda, tensor pres, tensor virial,
                   real *prescorr, real *enercorr, real *dvdlcorr);

void initialize_lambdas(FILE *fplog, t_inputrec *ir, int *fep_state, real *lambda, double *lam0);

void do_constrain_first(FILE *log, gmx_constr_t constr,
                        t_inputrec *inputrec, t_mdatoms *md,
                        t_state *state, t_commrec *cr, t_nrnb *nrnb,
                        t_forcerec *fr, gmx_localtop_t *top);

void init_md(FILE *fplog,
             t_commrec *cr, t_inputrec *ir, const output_env_t oenv,
             double *t, double *t0,
             real *lambda, int *fep_state, double *lam0,
             t_nrnb *nrnb, gmx_mtop_t *mtop,
             gmx_update_t *upd,
             int nfile, const t_filenm fnm[],
             gmx_mdoutf_t *outf, t_mdebin **mdebin,
             tensor force_vir, tensor shake_vir,
             rvec mu_tot,
             gmx_bool *bSimAnn, t_vcm **vcm, unsigned long Flags,
             gmx_wallcycle_t wcycle);
/* Routine in sim_util.c */

gmx_bool use_GPU(const struct nonbonded_verlet_t *nbv);

#ifdef __cplusplus
}
#endif

#endif  /* _sim_util_h */
