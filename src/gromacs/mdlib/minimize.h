/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "config.h"

#include <math.h>
#include <string.h>
#include <time.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/mtxio.h"
#include "gromacs/fileio/trajectory_writing.h"
#include "gromacs/imd/imd.h"
#include "gromacs/legacyheaders/constr.h"
#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/md_logging.h"
#include "gromacs/legacyheaders/md_support.h"
#include "gromacs/legacyheaders/mdatoms.h"
#include "gromacs/legacyheaders/mdebin.h"
#include "gromacs/legacyheaders/mdrun.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/legacyheaders/ns.h"
#include "gromacs/legacyheaders/sim_util.h"
#include "gromacs/legacyheaders/tgroup.h"
#include "gromacs/legacyheaders/txtdump.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/update.h"
#include "gromacs/legacyheaders/vsite.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/linearalgebra/sparsematrix.h"
#include "gromacs/listed-forces/manage-threading.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

typedef struct {
    t_state  s;
    rvec    *f;
    real     epot;
    real     fnorm;
    real     fmax;
    int      a_fmax;
} em_state_t;

em_state_t *init_em_state();

void print_em_start(FILE                     *fplog,
                    t_commrec                *cr,
                    gmx_walltime_accounting_t walltime_accounting,
                    gmx_wallcycle_t           wcycle,
                    const char               *name);

void em_time_end(gmx_walltime_accounting_t walltime_accounting,
                 gmx_wallcycle_t           wcycle);

void sp_header(FILE *out, const char *minimizer, real ftol, int nsteps);

void warn_step(FILE *fp, real ftol, gmx_bool bLastStep, gmx_bool bConstrain);

void print_converged(FILE *fp, const char *alg, real ftol,
                     gmx_int64_t count, gmx_bool bDone, gmx_int64_t nsteps,
                     real epot, real fmax, int nfmax, real fnorm);

void get_f_norm_max(t_commrec *cr,
                    t_grpopts *opts, t_mdatoms *mdatoms, rvec *f,
                    real *fnorm, real *fmax, int *a_fmax);

void get_state_f_norm_max(t_commrec *cr,
                          t_grpopts *opts, t_mdatoms *mdatoms,
                          em_state_t *ems);

void init_em(FILE *fplog, const char *title,
             t_commrec *cr, t_inputrec *ir,
             t_state *state_global, gmx_mtop_t *top_global,
             em_state_t *ems, gmx_localtop_t **top,
             rvec **f, rvec **f_global,
             t_nrnb *nrnb, rvec mu_tot,
             t_forcerec *fr, gmx_enerdata_t **enerd,
             t_graph **graph, t_mdatoms *mdatoms, gmx_global_stat_t *gstat,
             gmx_vsite_t *vsite, gmx_constr_t constr,
             int nfile, const t_filenm fnm[],
             gmx_mdoutf_t *outf, t_mdebin **mdebin,
             int imdport, unsigned long gmx_unused Flags,
             gmx_wallcycle_t wcycle);

void finish_em(t_commrec *cr, gmx_mdoutf_t outf,
               gmx_walltime_accounting_t walltime_accounting,
               gmx_wallcycle_t wcycle);

void swap_em_state(em_state_t *ems1, em_state_t *ems2);

void copy_em_coords(em_state_t *ems, t_state *state);

void write_em_traj(FILE *fplog, t_commrec *cr,
                   gmx_mdoutf_t outf,
                   gmx_bool bX, gmx_bool bF, const char *confout,
                   gmx_mtop_t *top_global,
                   t_inputrec *ir, gmx_int64_t step,
                   em_state_t *state,
                   t_state *state_global, rvec *f_global);

void do_em_step(t_commrec *cr, t_inputrec *ir, t_mdatoms *md,
                gmx_bool bMolPBC,
                em_state_t *ems1, real a, rvec *f, em_state_t *ems2,
                gmx_constr_t constr, gmx_localtop_t *top,
                t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                gmx_int64_t count);

void em_dd_partition_system(FILE *fplog, int step, t_commrec *cr,
                            gmx_mtop_t *top_global, t_inputrec *ir,
                            em_state_t *ems, gmx_localtop_t *top,
                            t_mdatoms *mdatoms, t_forcerec *fr,
                            gmx_vsite_t *vsite, gmx_constr_t constr,
                            t_nrnb *nrnb, gmx_wallcycle_t wcycle);

void evaluate_energy(FILE *fplog, t_commrec *cr,
                     gmx_mtop_t *top_global,
                     em_state_t *ems, gmx_localtop_t *top,
                     t_inputrec *inputrec,
                     t_nrnb *nrnb, gmx_wallcycle_t wcycle,
                     gmx_global_stat_t gstat,
                     gmx_vsite_t *vsite, gmx_constr_t constr,
                     t_fcdata *fcd,
                     t_graph *graph, t_mdatoms *mdatoms,
                     t_forcerec *fr, rvec mu_tot,
                     gmx_enerdata_t *enerd, tensor vir, tensor pres,
                     gmx_int64_t count, gmx_bool bFirst);

double reorder_partsum(t_commrec *cr, t_grpopts *opts, t_mdatoms *mdatoms,
                       gmx_mtop_t *mtop,
                       em_state_t *s_min, em_state_t *s_b);

real pr_beta(t_commrec *cr, t_grpopts *opts, t_mdatoms *mdatoms,
             gmx_mtop_t *mtop,
             em_state_t *s_min, em_state_t *s_b);

} /* namespace gmx */
