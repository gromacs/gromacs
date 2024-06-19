/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_MDLIB_SHELLFC_H
#define GMX_MDLIB_SHELLFC_H

#include <cstdint>
#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

class DDBalanceRegionHandler;
struct gmx_enerdata_t;
struct gmx_enfrot;
struct gmx_localtop_t;
struct gmx_multisim_t;
struct gmx_shellfc_t;
struct gmx_mtop_t;
class history_t;
struct pull_t;
struct t_forcerec;
struct t_inputrec;
struct t_mdatoms;
struct t_nrnb;
class t_state;
class CpuPpLongRangeNonbondeds;
struct t_commrec;

namespace gmx
{
template<typename>
class ArrayRef;
template<typename>
class ArrayRefWithPadding;
class Constraints;
class ForceBuffersView;
class ImdSession;
struct MDModulesNotifiers;
class MdrunScheduleWorkload;
class VirtualSitesHandler;
} // namespace gmx

/*! \brief Initialization function, also predicts the initial shell positions.
 *
 * \param fplog Pointer to the log stream. Can be set to \c nullptr to disable verbose log.
 * \param mtop Pointer to a global system topology object.
 * \param nflexcon Number of flexible constraints.
 * \param nstcalcenergy How often are energies calculated. Must be provided for sanity check.
 * \param usingDomainDecomposition Whether domain decomposition is used. Must be provided for sanity check.
 * \param usingPmeOnGpu Set to true if GPU will be used for PME calculations. Necessary for proper buffer initialization.
 *
 * \returns a pointer to an initialized \c shellfc object.
 */
gmx_shellfc_t* init_shell_flexcon(FILE*             fplog,
                                  const gmx_mtop_t& mtop,
                                  int               nflexcon,
                                  int               nstcalcenergy,
                                  bool              usingDomainDecomposition,
                                  bool              usingPmeOnGpu);

/* Optimize shell positions */
void relax_shell_flexcon(FILE*                               log,
                         const t_commrec*                    cr,
                         const gmx_multisim_t*               ms,
                         gmx_bool                            bVerbose,
                         gmx_enfrot*                         enforcedRotation,
                         int64_t                             mdstep,
                         const t_inputrec*                   inputrec,
                         const gmx::MDModulesNotifiers&      mdModulesNotifiers,
                         gmx::ImdSession*                    imdSession,
                         pull_t*                             pull_work,
                         gmx_bool                            bDoNS,
                         const gmx_localtop_t*               top,
                         gmx::Constraints*                   constr,
                         gmx_enerdata_t*                     enerd,
                         int                                 natoms,
                         gmx::ArrayRefWithPadding<gmx::RVec> x,
                         gmx::ArrayRefWithPadding<gmx::RVec> v,
                         const matrix                        box,
                         gmx::ArrayRef<real>                 lambda,
                         const history_t*                    hist,
                         gmx::ForceBuffersView*              f,
                         tensor                              force_vir,
                         const t_mdatoms&                    md,
                         CpuPpLongRangeNonbondeds*           longRangeNonbondeds,
                         t_nrnb*                             nrnb,
                         gmx_wallcycle*                      wcycle,
                         gmx_shellfc_t*                      shfc,
                         t_forcerec*                         fr,
                         const gmx::MdrunScheduleWorkload&   runScheduleWork,
                         double                              t,
                         rvec                                mu_tot,
                         gmx::VirtualSitesHandler*           vsite,
                         const DDBalanceRegionHandler&       ddBalanceRegionHandler);

/* Print some final output and delete shellfc */
void done_shellfc(FILE* fplog, gmx_shellfc_t* shellfc, int64_t numSteps);

#endif
