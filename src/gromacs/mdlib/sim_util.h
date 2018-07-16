/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_SIM_UTIL_H
#define GMX_MDLIB_SIM_UTIL_H

#include <gromacs/mdtypes/iforceprovider.h>
#include "gromacs/fileio/enxio.h"
#include "gromacs/mdlib/mdebin.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gmx_omp_nthreads.h"

struct ForceProviders;
struct gmx_output_env_t;
struct gmx_pme_t;
struct gmx_update_t;
struct MdrunOptions;
struct nonbonded_verlet_t;
struct t_forcerec;
struct t_mdatoms;
struct t_nrnb;
struct gmx_vsite_t;
struct interaction_const_t;

struct gmx_enfrot;
struct gmx_edsam;

namespace gmx
{
class BoxDeformation;
class Constraints;
class IMDOutputProvider;
class MDLogger;
class ForceWithVirial;
}

typedef struct gmx_global_stat *gmx_global_stat_t;

void do_pbc_first_mtop(FILE *fplog, int ePBC, const matrix box,
                       const gmx_mtop_t *mtop, rvec x[]);

void do_pbc_mtop(FILE *fplog, int ePBC, const matrix box,
                 const gmx_mtop_t *mtop, rvec x[]);

/*! \brief Parallellizes put_atoms_in_box()
 *
 * This wrapper function around put_atoms_in_box() with the ugly manual
 * workload splitting is needed to avoid silently introducing multithreading
 * in tools.
 * \param[in]    ePBC   The pbc type
 * \param[in]    box    The simulation box
 * \param[inout] x      The coordinates of the atoms
 */
void put_atoms_in_box_omp(int ePBC, const matrix box, gmx::ArrayRef<gmx::RVec> x);



/* ROUTINES from stat.c */
gmx_global_stat_t global_stat_init(const t_inputrec *ir);

void global_stat_destroy(gmx_global_stat_t gs);

void global_stat(const gmx_global_stat *gs,
                 const t_commrec *cr, gmx_enerdata_t *enerd,
                 tensor fvir, tensor svir, rvec mu_tot,
                 const t_inputrec *inputrec,
                 gmx_ekindata_t *ekind,
                 const gmx::Constraints *constr, t_vcm *vcm,
                 int nsig, real *sig,
                 int *totalNumberOfBondedInteractions,
                 gmx_bool bSumEkinhOld, int flags);
/* All-reduce energy-like quantities over cr->mpi_comm_mysim */

bool do_per_step(int64_t step, int64_t nstep);
/* Return TRUE if io should be done */

/* ROUTINES from sim_util.c */

void print_time(FILE *out, gmx_walltime_accounting_t walltime_accounting,
                int64_t step, const t_inputrec *ir, const t_commrec *cr);

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

void print_start(FILE *fplog, const t_commrec *cr,
                 gmx_walltime_accounting_t walltime_accounting,
                 const char *name);

void finish_run(FILE *log, const gmx::MDLogger &mdlog, const t_commrec *cr,
                const t_inputrec *inputrec,
                t_nrnb nrnb[], gmx_wallcycle_t wcycle,
                gmx_walltime_accounting_t walltime_accounting,
                nonbonded_verlet_t *nbv,
                const gmx_pme_t *pme,
                gmx_bool bWriteStat);

void calc_enervirdiff(FILE *fplog, int eDispCorr, t_forcerec *fr);

void calc_dispcorr(const t_inputrec *ir, const t_forcerec *fr,
                   const matrix box, real lambda, tensor pres, tensor virial,
                   real *prescorr, real *enercorr, real *dvdlcorr);

void initialize_lambdas(FILE *fplog, t_inputrec *ir, int *fep_state, gmx::ArrayRef<real> lambda, double *lam0);

void do_constrain_first(FILE *log, gmx::Constraints *constr,
                        const t_inputrec *inputrec, const t_mdatoms *md,
                        t_state *state);

void init_md(FILE *fplog,
             const t_commrec *cr, gmx::IMDOutputProvider *outputProvider,
             t_inputrec *ir, const gmx_output_env_t *oenv,
             const MdrunOptions &mdrunOptions,
             double *t, double *t0,
             t_state *globalState, double *lam0,
             t_nrnb *nrnb, gmx_mtop_t *mtop,
             gmx_update_t **upd,
             gmx::BoxDeformation *deform,
             int nfile, const t_filenm fnm[],
             gmx_mdoutf_t *outf, t_mdebin **mdebin,
             tensor force_vir, tensor shake_vir,
             tensor total_vir, tensor pres,
             rvec mu_tot,
             gmx_bool *bSimAnn, t_vcm **vcm,
             gmx_wallcycle_t wcycle);

void init_rerun(FILE *fplog,
                const t_commrec *cr, gmx::IMDOutputProvider *outputProvider,
                t_inputrec *ir, const gmx_output_env_t *oenv,
                const MdrunOptions &mdrunOptions,
                t_state *globalState, double *lam0,
                t_nrnb *nrnb, gmx_mtop_t *mtop,
                int nfile, const t_filenm fnm[],
                gmx_mdoutf_t *outf, t_mdebin **mdebin,
                gmx_wallcycle_t wcycle);

/* Routine in sim_util.c */

gmx_bool use_GPU(const nonbonded_verlet_t *nbv);

void post_process_forces(const t_commrec           *cr,
                         int64_t                    step,
                         t_nrnb                    *nrnb,
                         gmx_wallcycle_t            wcycle,
                         const gmx_localtop_t      *top,
                         const matrix               box,
                         const rvec                 x[],
                         rvec                       f[],
                         gmx::ForceWithVirial      *forceWithVirial,
                         tensor                     vir_force,
                         const t_mdatoms           *mdatoms,
                         const t_graph             *graph,
                         const t_forcerec          *fr,
                         const gmx_vsite_t         *vsite,
                         int                        flags);

inline void clear_rvecs_omp(int n, rvec v[])
{
    int nth = gmx_omp_nthreads_get_simple_rvec_task(emntDefault, n);

    /* Note that we would like to avoid this conditional by putting it
     * into the omp pragma instead, but then we still take the full
     * omp parallel for overhead (at least with gcc5).
     */
    if (nth == 1)
    {
        for (int i = 0; i < n; i++)
        {
            clear_rvec(v[i]);
        }
    }
    else
    {
#pragma omp parallel for num_threads(nth) schedule(static)
        for (int i = 0; i < n; i++)
        {
            clear_rvec(v[i]);
        }
    }
}

void do_nb_verlet(const t_forcerec *fr,
                  const interaction_const_t *ic,
                  gmx_enerdata_t *enerd,
                  int flags, int ilocality,
                  int clearF,
                  int64_t step,
                  t_nrnb *nrnb,
                  gmx_wallcycle_t wcycle);

void do_nb_verlet_fep(nbnxn_pairlist_set_t *nbl_lists,
                      t_forcerec           *fr,
                      rvec                  x[],
                      rvec                  f[],
                      const t_mdatoms      *mdatoms,
                      t_lambda             *fepvals,
                      real                 *lambda,
                      gmx_enerdata_t       *enerd,
                      int                   flags,
                      t_nrnb               *nrnb,
                      gmx_wallcycle_t       wcycle);

void calc_virial(int start, int homenr, const rvec x[], const rvec f[],
                 tensor vir_part, const t_graph *graph, const matrix box,
                 t_nrnb *nrnb, const t_forcerec *fr, int ePBC);

void checkPotentialEnergyValidity(int64_t               step,
                                  const gmx_enerdata_t &enerd,
                                  const t_inputrec     &inputrec);

void
computeSpecialForces(FILE                          *fplog,
                     const t_commrec               *cr,
                     const t_inputrec              *inputrec,
                     gmx::Awh                      *awh,
                     gmx_enfrot                    *enforcedRotation,
                     int64_t                        step,
                     double                         t,
                     gmx_wallcycle_t                wcycle,
                     ForceProviders                *forceProviders,
                     matrix                         box,
                     gmx::ArrayRef<const gmx::RVec> x,
                     const t_mdatoms               *mdatoms,
                     real                          *lambda,
                     int                            forceFlags,
                     gmx::ForceWithVirial          *forceWithVirial,
                     gmx_enerdata_t                *enerd,
                     gmx_edsam                     *ed,
                     gmx_bool                       bNS);

#endif
