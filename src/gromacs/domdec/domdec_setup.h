/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*! \internal \file
 *
 * \brief Declares functions for choosing the DD grid setup
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */
#ifndef GMX_DOMDEC_DOMDEC_SETUP_H
#define GMX_DOMDEC_DOMDEC_SETUP_H


#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/real.h"

enum class DDRole;
struct DDSettings;
struct DDSystemInfo;
struct gmx_ddbox_t;
struct gmx_mtop_t;
struct t_commrec;
struct t_inputrec;

namespace gmx
{
struct DomdecOptions;
struct MDModulesNotifiers;
class MDLogger;
class SeparatePmeRanksPermitted;
template<typename T>
class ArrayRef;
} // namespace gmx

/*! \brief Returns the volume fraction of the system that is communicated */
real comm_box_frac(const gmx::IVec& dd_nc, real cutoff, const gmx_ddbox_t& ddbox);

/*! \internal
 * \brief Describes the DD grid setup
 *
 * This struct is for temporary use when choosing and initializing
 * the domain decomposition grid.
 */
struct DDGridSetup
{
    //! The number of separate PME ranks, 0 if none or all ranks do PME
    int numPmeOnlyRanks = 0;
    //! The number of domains along each dimension
    ivec numDomains = { 0, 0, 0 };
    //! The number of dimensions which we decompose in domains
    int numDDDimensions = 0;
    //! The domain decomposition dimensions, the first numDDDimensions entries are used
    ivec ddDimensions = { -1, -1, -1 };
};

/*! \brief Checks for ability to use separate PME ranks
 *
 * Disables automatic usage if:
 * some MDModule could not use separate PME ranks,
 * GPU setup is not compatible with separate PME ranks,
 * user provided explicit DD grid
 * or total number of ranks is not large enough to use PME ranks
 */
gmx::SeparatePmeRanksPermitted checkForSeparatePmeRanks(const gmx::MDModulesNotifiers& notifiers,
                                                        const gmx::DomdecOptions&      options,
                                                        int  numRanksRequested,
                                                        bool useGpuForNonbonded,
                                                        bool useGpuForPme,
                                                        bool canUseGpuPmeDecomposition);

/*! \brief Checks that requests for PP and PME ranks honor basic expectations
 *
 * Issues a fatal error if there are more PME ranks than PP, if the
 * count of PP ranks has a prime factor that is too large to be likely
 * to have good performance or PME-only ranks could not be used,
 * but requested with -npme > 0 */
void checkForValidRankCountRequests(int                                   numRanksRequested,
                                    bool                                  usingPme,
                                    int                                   numPmeRanksRequested,
                                    const gmx::SeparatePmeRanksPermitted& separatePmeRanksPermitted,
                                    bool checkForLargePrimeFactors);

/*! \brief Return the minimum cell size (in nm) required for DD */
real getDDGridSetupCellSizeLimit(const gmx::MDLogger& mdlog,
                                 bool                 bDynLoadBal,
                                 real                 dlb_scale,
                                 const t_inputrec&    ir,
                                 real                 systemInfoCellSizeLimit,
                                 int                  numRanksRequested);

/*! \brief Determines the DD grid setup
 *
 * Either implements settings required by the user, or otherwise
 * chooses estimated optimal number of separate PME ranks and DD grid
 * cell setup, DD cell size limits, and the initial ddbox.
 */
DDGridSetup getDDGridSetup(const gmx::MDLogger&                  mdlog,
                           DDRole                                ddRole,
                           MPI_Comm                              communicator,
                           int                                   numRanksRequested,
                           const gmx::DomdecOptions&             options,
                           const DDSettings&                     ddSettings,
                           const DDSystemInfo&                   systemInfo,
                           real                                  cellSizeLimit,
                           const gmx_mtop_t&                     mtop,
                           const t_inputrec&                     ir,
                           const gmx::SeparatePmeRanksPermitted& separatePmeRanksPermitted,
                           const matrix                          box,
                           gmx::ArrayRef<const gmx::RVec>        xGlobal,
                           gmx_ddbox_t*                          ddbox);

#endif
