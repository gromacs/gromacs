/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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

#ifndef GMX_MDTYPES_DEFAULTSCHEDULE_H
#define GMX_MDTYPES_DEFAULTSCHEDULE_H

#include "gromacs/math/vectypes.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/timing/wallcycle.h"

struct t_commrec;
struct t_inputrec;
struct t_forcerec;
struct t_mdatoms;
struct t_nrnb;
struct gmx_edsam;
struct gmx_enerdata_t;
struct gmx_enfrot;
struct gmx_groups_t;
struct gmx_localtop_t;
struct gmx_multisim_t;
struct gmx_vsite_t;
struct t_fcdata;
struct t_graph;
struct interaction_const_t;

class history_t;

enum class DdOpenBalanceRegionBeforeForceComputation;
enum class DdCloseBalanceRegionAfterForceComputation;

namespace gmx
{
class PpForceWorkload;
class Awh;
}   // namespace gmx

void do_force_cutsVERLET(FILE *fplog,
                         const t_commrec *cr,
                         const gmx_multisim_t *ms,
                         const t_inputrec *inputrec,
                         gmx::Awh *awh,
                         gmx_enfrot *enforcedRotation,
                         int64_t step,
                         t_nrnb *nrnb,
                         gmx_wallcycle_t wcycle,
                         const gmx_localtop_t *top,
                         const gmx_groups_t * /* groups */,
                         matrix box, gmx::ArrayRefWithPadding<gmx::RVec> x,
                         history_t *hist,
                         gmx::ArrayRefWithPadding<gmx::RVec> force,
                         tensor vir_force,
                         const t_mdatoms *mdatoms,
                         gmx_enerdata_t *enerd, t_fcdata *fcd,
                         real *lambda,
                         t_graph *graph,
                         t_forcerec *fr,
                         gmx::PpForceWorkload *ppForceWorkload,
                         interaction_const_t *ic,
                         const gmx_vsite_t *vsite,
                         rvec mu_tot,
                         double t,
                         gmx_edsam *ed,
                         int flags,
                         DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion,
                         DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion);

void do_force_cutsGROUP(FILE *fplog,
                        const t_commrec *cr,
                        const gmx_multisim_t *ms,
                        const t_inputrec *inputrec,
                        gmx::Awh *awh,
                        gmx_enfrot *enforcedRotation,
                        int64_t step,
                        t_nrnb *nrnb,
                        gmx_wallcycle_t wcycle,
                        gmx_localtop_t *top,
                        const gmx_groups_t *groups,
                        matrix box, gmx::ArrayRefWithPadding<gmx::RVec> x,
                        history_t *hist,
                        gmx::ArrayRefWithPadding<gmx::RVec> force,
                        tensor vir_force,
                        const t_mdatoms *mdatoms,
                        gmx_enerdata_t *enerd,
                        t_fcdata *fcd,
                        real *lambda,
                        t_graph *graph,
                        t_forcerec *fr,
                        const gmx_vsite_t *vsite,
                        rvec mu_tot,
                        double t,
                        gmx_edsam *ed,
                        int flags,
                        DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion,
                        DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion);

void do_force(FILE                                     *fplog,
              const t_commrec                          *cr,
              const gmx_multisim_t                     *ms,
              const t_inputrec                         *inputrec,
              gmx::Awh                                 *awh,
              gmx_enfrot                               *enforcedRotation,
              int64_t                                   step,
              t_nrnb                                   *nrnb,
              gmx_wallcycle_t                           wcycle,
              gmx_localtop_t                           *top,
              const gmx_groups_t                       *groups,
              matrix                                    box,
              gmx::ArrayRefWithPadding<gmx::RVec>       x,     //NOLINT(performance-unnecessary-value-param)
              history_t                                *hist,
              gmx::ArrayRefWithPadding<gmx::RVec>       force, //NOLINT(performance-unnecessary-value-param)
              tensor                                    vir_force,
              const t_mdatoms                          *mdatoms,
              gmx_enerdata_t                           *enerd,
              t_fcdata                                 *fcd,
              gmx::ArrayRef<real>                       lambda,
              t_graph                                  *graph,
              t_forcerec                               *fr,
              gmx::PpForceWorkload                     *ppForceWorkload,
              const gmx_vsite_t                        *vsite,
              rvec                                      mu_tot,
              double                                    t,
              gmx_edsam                                *ed,
              int                                       flags,
              DdOpenBalanceRegionBeforeForceComputation ddOpenBalanceRegion,
              DdCloseBalanceRegionAfterForceComputation ddCloseBalanceRegion);

#endif //GMX_MDTYPES_DEFAULTSCHEDULE_H
