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
#ifndef GMX_MDLIB_FORCE_H
#define GMX_MDLIB_FORCE_H

#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"

struct gmx_edsam;
struct gmx_enerdata_t;
struct gmx_groups_t;
struct gmx_grppairener_t;
struct gmx_localtop_t;
struct gmx_multisim_t;
struct gmx_vsite_t;
struct gmx_wallcycle;
class history_t;
struct t_blocka;
struct t_commrec;
struct t_fcdata;
struct t_forcerec;
struct t_graph;
struct t_idef;
struct t_inputrec;
struct t_lambda;
struct t_mdatoms;
struct t_nrnb;

namespace gmx
{
class Awh;
class ForceWithVirial;
class MDLogger;
}

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

void do_force(FILE                                     *log,
              const t_commrec                          *cr,
              const gmx_multisim_t                     *ms,
              const t_inputrec                         *inputrec,
              gmx::Awh                                 *awh,
              int64_t                                   step,
              t_nrnb                                   *nrnb,
              gmx_wallcycle                            *wcycle,
              // TODO top can be const when the group scheme no longer
              // builds exclusions during neighbor searching within
              // do_force_cutsGROUP.
              gmx_localtop_t                           *top,
              const gmx_groups_t                       *groups,
              matrix                                    box,
              gmx::PaddedArrayRef<gmx::RVec>            coordinates,
              history_t                                *hist,
              gmx::PaddedArrayRef<gmx::RVec>            force,
              tensor                                    vir_force,
              const t_mdatoms                          *mdatoms,
              gmx_enerdata_t                           *enerd,
              t_fcdata                                 *fcd,
              gmx::ArrayRef<real>                       lambda,
              t_graph                                  *graph,
              t_forcerec                               *fr,
              const gmx_vsite_t                        *vsite,
              rvec                                      mu_tot,
              double                                    t,
              const gmx_edsam                          *ed,
              int                                       flags,
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

void ns(FILE               *fplog,
        t_forcerec         *fr,
        matrix              box,
        const gmx_groups_t *groups,
        gmx_localtop_t     *top,
        const t_mdatoms    *md,
        const t_commrec    *cr,
        t_nrnb             *nrnb,
        gmx_bool            bFillGrid);
/* Call the neighborsearcher */

void do_force_lowlevel(t_forcerec   *fr,
                       const t_inputrec *ir,
                       const t_idef *idef,
                       const t_commrec *cr,
                       const gmx_multisim_t *ms,
                       t_nrnb       *nrnb,
                       gmx_wallcycle *wcycle,
                       const t_mdatoms *md,
                       rvec         x[],
                       history_t    *hist,
                       rvec         f_shortrange[],
                       gmx::ForceWithVirial *forceWithVirial,
                       gmx_enerdata_t *enerd,
                       t_fcdata     *fcd,
                       matrix       box,
                       t_lambda     *fepvals,
                       real         *lambda,
                       const t_graph *graph,
                       const t_blocka *excl,
                       rvec         mu_tot[2],
                       int          flags,
                       float        *cycles_pme);
/* Call all the force routines */

#endif
