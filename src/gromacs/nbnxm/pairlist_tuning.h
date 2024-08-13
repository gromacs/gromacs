/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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

/*! \libinternal \file
 *
 * \brief Declares functions for tuning adjustable parameters for the nbnxn non-bonded search and interaction kernels
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \inlibraryapi
 * \ingroup module_nbnxm
 */

#ifndef NBNXM_PAIRLIST_TUNING_H
#define NBNXM_PAIRLIST_TUNING_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct interaction_const_t;
struct t_commrec;
struct t_inputrec;

namespace gmx
{
struct PairlistParams;
class CpuInfo;
class MDLogger;

/*! \brief Try to increase nstlist when using the Verlet cut-off scheme
 *
 * \param[in,out] fplog    Log file
 * \param[in]     cr       The communication record
 * \param[in]     ir       The input parameter record
 * \param[in]     nstlistOnCmdline  The value of nstlist provided on the command line
 * \param[in]     mtop     The global topology
 * \param[in]     box      The unit cell
 * \param[in]     effectiveAtomDensity  The effective atom density
 * \param[in]     useOrEmulateGpuForNonbondeds  Tells if we are using a GPU for non-bondeds
 * \param[in]     cpuinfo  Information about the CPU(s)
 */
void increaseNstlist(FILE*             fplog,
                     t_commrec*        cr,
                     t_inputrec*       ir,
                     int               nstlistOnCmdline,
                     const gmx_mtop_t* mtop,
                     const matrix      box,
                     real              effectiveAtomDensity,
                     bool              useOrEmulateGpuForNonbondeds,
                     const CpuInfo&    cpuinfo);

/*! \brief Set up the dynamic pairlist pruning
 *
 * \param[in,out] mdlog            MD logger
 * \param[in]     inputrec         The input parameter record
 * \param[in]     mtop             The global topology
 * \param[in]     effectiveAtomDensity  The effective atom density of the system
 * \param[in]     interactionConst The nonbonded interactions constants
 * \param[in,out] listParams       The list setup parameters
 */
void setupDynamicPairlistPruning(const MDLogger&            mdlog,
                                 const t_inputrec&          inputrec,
                                 const gmx_mtop_t&          mtop,
                                 real                       effectiveAtomDensity,
                                 const interaction_const_t& interactionConst,
                                 PairlistParams*            listParams);

/*! \brief Prints an estimate of the error in the pressure due to missing interactions
 *
 * The NBNxM algorithm tolerates a few missing pair interactions.
 * Missing pair interactions will lead to a systematic overestimates of
 * the pressure when dispersion forces dominate at the cut-off distance.
 * This routine prints an overestimate of the error in the average pressure.
 *
 * \param[in,out] mdlog            MD logger
 * \param[in]     inputrec         The input parameter record
 * \param[in]     mtop             The global topology
 * \param[in]     effectiveAtomDensity  The effective atom density of the system
 * \param[in]     listParams       The list setup parameters
 */
void printNbnxmPressureError(const MDLogger&       mdlog,
                             const t_inputrec&     inputrec,
                             const gmx_mtop_t&     mtop,
                             real                  effectiveAtomDensity,
                             const PairlistParams& listParams);

} // namespace gmx

#endif /* NBNXM_PAIRLIST_TUNING_H */
