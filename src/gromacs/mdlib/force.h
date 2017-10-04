/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_FORCE_H
#define GMX_MDLIB_FORCE_H

#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arrayref.h"

struct gmx_device_info_t;
struct gmx_edsam;
struct gmx_gpu_info_t;
struct gmx_groups_t;
struct gmx_vsite_t;
class history_t;
struct nonbonded_verlet_t;
struct t_blocka;
struct t_commrec;
struct t_fcdata;
struct t_graph;
struct t_grpopts;
struct t_inputrec;
struct t_lambda;
struct t_mdatoms;
struct t_nrnb;
struct t_pbc;

namespace gmx
{
class ForceWithVirial;
class MDLogger;
}

void calc_vir(int nxf, rvec x[], rvec f[], tensor vir,
              gmx_bool bScrewPBC, matrix box);
/* Calculate virial for nxf atoms, and add it to vir */

void f_calc_vir(int i0, int i1, rvec x[], rvec f[], tensor vir,
                t_graph *g, rvec shift_vec[]);
/* Calculate virial taking periodicity into account */

real RF_excl_correction(const t_forcerec *fr, t_graph *g,
                        const t_mdatoms *mdatoms, const t_blocka *excl,
                        rvec x[], rvec f[], rvec *fshift, const t_pbc *pbc,
                        real lambda, real *dvdlambda);
/* Calculate the reaction-field energy correction for this node:
 * epsfac q_i q_j (k_rf r_ij^2 - c_rf)
 * and force correction for all excluded pairs, including self pairs.
 */

void calc_rffac(FILE *fplog, int eel, real eps_r, real eps_rf,
                real Rc, real Temp,
                real zsq, matrix box,
                real *krf, real *crf);
/* Determine the reaction-field constants */

void init_generalized_rf(FILE *fplog,
                         const gmx_mtop_t *mtop, const t_inputrec *ir,
                         t_forcerec *fr);
/* Initialize the generalized reaction field parameters */


/* In wall.c */
void make_wall_tables(FILE *fplog,
                      const t_inputrec *ir, const char *tabfn,
                      const gmx_groups_t *groups,
                      t_forcerec *fr);

real do_walls(t_inputrec *ir, t_forcerec *fr, matrix box, t_mdatoms *md,
              rvec x[], rvec f[], real lambda, real Vlj[], t_nrnb *nrnb);

gmx_bool can_use_allvsall(const t_inputrec *ir,
                          gmx_bool bPrintNote, t_commrec *cr, FILE *fp);
/* Returns if we can use all-vs-all loops.
 * If bPrintNote==TRUE, prints a note, if necessary, to stderr
 * and fp (if !=NULL) on the master node.
 */

gmx_bool nbnxn_simd_supported(const gmx::MDLogger &mdlog,
                              const t_inputrec    *ir);
/* Return if CPU SIMD support exists for the given inputrec
 * If the return value is FALSE and fplog/cr != NULL, prints a fallback
 * message to fplog/stderr.
 */

gmx_bool uses_simple_tables(int                 cutoff_scheme,
                            nonbonded_verlet_t *nbv,
                            int                 group);
/* Returns whether simple tables (i.e. not for use with GPUs) are used
 * with the type of kernel indicated.
 */

void init_enerdata(int ngener, int n_lambda, gmx_enerdata_t *enerd);
/* Intializes the energy storage struct */

void destroy_enerdata(gmx_enerdata_t *enerd);
/* Free all memory associated with enerd */

void reset_foreign_enerdata(gmx_enerdata_t *enerd);
/* Resets only the foreign energy data */

void reset_enerdata(gmx_enerdata_t *enerd);
/* Resets the energy data */

void sum_epot(gmx_grppairener_t *grpp, real *epot);
/* Locally sum the non-bonded potential energy terms */

void sum_dhdl(gmx_enerdata_t *enerd, gmx::ArrayRef<const real> lambda, t_lambda *fepvals);
/* Sum the free energy contributions */

/* Compute the average C6 and C12 params for LJ corrections */
void set_avcsixtwelve(FILE *fplog, t_forcerec *fr,
                      const gmx_mtop_t *mtop);

void do_force(FILE *log, t_commrec *cr,
              t_inputrec *inputrec,
              gmx_int64_t step, struct t_nrnb *nrnb, gmx_wallcycle_t wcycle,
              gmx_localtop_t *top,
              gmx_groups_t *groups,
              matrix box, PaddedRVecVector *coordinates, history_t *hist,
              PaddedRVecVector *force,
              tensor vir_force,
              t_mdatoms *mdatoms,
              gmx_enerdata_t *enerd, t_fcdata *fcd,
              gmx::ArrayRef<real> lambda, t_graph *graph,
              t_forcerec *fr,
              gmx_vsite_t *vsite, rvec mu_tot,
              double t, struct gmx_edsam *ed,
              gmx_bool bBornRadii,
              int flags,
              DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion,
              DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion);

/* Communicate coordinates (if parallel).
 * Do neighbor searching (if necessary).
 * Calculate forces.
 * Communicate forces (if parallel).
 * Spread forces for vsites (if present).
 *
 * f is always required.
 */

void ns(FILE              *fplog,
        t_forcerec        *fr,
        matrix             box,
        gmx_groups_t      *groups,
        gmx_localtop_t    *top,
        t_mdatoms         *md,
        t_commrec         *cr,
        t_nrnb            *nrnb,
        gmx_bool           bFillGrid);
/* Call the neighborsearcher */

void do_force_lowlevel(t_forcerec   *fr,
                       t_inputrec   *ir,
                       t_idef       *idef,
                       t_commrec    *cr,
                       t_nrnb       *nrnb,
                       gmx_wallcycle_t wcycle,
                       t_mdatoms    *md,
                       rvec         x[],
                       history_t    *hist,
                       rvec         f_shortrange[],
                       gmx::ForceWithVirial *forceWithVirial,
                       gmx_enerdata_t *enerd,
                       t_fcdata     *fcd,
                       gmx_localtop_t *top,
                       gmx_genborn_t *born,
                       gmx_bool         bBornRadii,
                       matrix       box,
                       t_lambda     *fepvals,
                       real         *lambda,
                       t_graph      *graph,
                       t_blocka     *excl,
                       rvec         mu_tot[2],
                       int          flags,
                       float        *cycles_pme);
/* Call all the force routines */

void free_gpu_resources(const t_forcerec            *fr,
                        const t_commrec             *cr,
                        const gmx_device_info_t     *deviceInfo);

#endif
