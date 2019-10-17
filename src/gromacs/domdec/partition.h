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
/*! \libinternal \file
 *
 * \brief This file declares functions for mdrun to call to make a new
 * domain decomposition, and check it.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_PARTITION_H
#define GMX_DOMDEC_PARTITION_H

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_ddbox_t;
struct gmx_domdec_t;
struct gmx_localtop_t;
struct gmx_mtop_t;
struct gmx_vsite_t;
struct gmx_wallcycle;
struct pull_t;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct t_nrnb;
class t_state;

namespace gmx
{
class Constraints;
class ImdSession;
class MDAtoms;
class MDLogger;
} // namespace gmx

//! Check whether the DD grid has moved too far for correctness.
bool check_grid_jump(int64_t step, const gmx_domdec_t* dd, real cutoff, const gmx_ddbox_t* ddbox, gmx_bool bFatal);

/*! \brief Print statistics for domain decomposition communication */
void print_dd_statistics(const t_commrec* cr, const t_inputrec* ir, FILE* fplog);

/*! \brief Partition the system over the nodes.
 *
 * step is only used for printing error messages.
 * If bMasterState==TRUE then state_global from the master node is used,
 * else state_local is redistributed between the nodes.
 * When f!=NULL, *f will be reallocated to the size of state_local.
 *
 * \param[in] fplog         Pointer to the log file
 * \param[in] mdlog         MD file logger
 * \param[in] step          Current step
 * \param[in] cr            Communication record
 * \param[in] bMasterState  Is it a master state
 * \param[in] nstglobalcomm Will globals be computed on this step
 * \param[in] state_global  Global state
 * \param[in] top_global    Global topology
 * \param[in] ir            Input record
 * \param[in] imdSession    IMD handle
 * \param[in] pull_work     Pulling data
 * \param[in] state_local   Local state
 * \param[in] f             Force buffer
 * \param[in] mdatoms       MD atoms
 * \param[in] top_local     Local topology
 * \param[in] fr            Force record
 * \param[in] vsite         Virtual sites
 * \param[in] constr        Constraints
 * \param[in] nrnb          Cycle counters
 * \param[in] wcycle        Timers
 * \param[in] bVerbose      Be verbose
 */
void dd_partition_system(FILE*                             fplog,
                         const gmx::MDLogger&              mdlog,
                         int64_t                           step,
                         const t_commrec*                  cr,
                         gmx_bool                          bMasterState,
                         int                               nstglobalcomm,
                         t_state*                          state_global,
                         const gmx_mtop_t&                 top_global,
                         const t_inputrec*                 ir,
                         gmx::ImdSession*                  imdSession,
                         pull_t*                           pull_work,
                         t_state*                          state_local,
                         gmx::PaddedHostVector<gmx::RVec>* f,
                         gmx::MDAtoms*                     mdatoms,
                         gmx_localtop_t*                   top_local,
                         t_forcerec*                       fr,
                         gmx_vsite_t*                      vsite,
                         gmx::Constraints*                 constr,
                         t_nrnb*                           nrnb,
                         gmx_wallcycle*                    wcycle,
                         gmx_bool                          bVerbose);

/*! \brief Check whether bonded interactions are missing, if appropriate
 *
 * \param[in]    mdlog                                  Logger
 * \param[in]    cr                                     Communication object
 * \param[in]    totalNumberOfBondedInteractions        Result of the global reduction over the number of bonds treated in each domain
 * \param[in]    top_global                             Global topology for the error message
 * \param[in]    top_local                              Local topology for the error message
 * \param[in]    x                                      Position vector for the error message
 * \param[in]    box                                    Box matrix for the error message
 * \param[in,out] shouldCheckNumberOfBondedInteractions Whether we should do the check. Always set to false.
 */
void checkNumberOfBondedInteractions(const gmx::MDLogger&  mdlog,
                                     t_commrec*            cr,
                                     int                   totalNumberOfBondedInteractions,
                                     const gmx_mtop_t*     top_global,
                                     const gmx_localtop_t* top_local,
                                     const rvec*           x,
                                     const matrix          box,
                                     bool*                 shouldCheckNumberOfBondedInteractions);

#endif
