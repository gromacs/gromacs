/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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

struct DDSettings;
struct DDSystemInfo;
struct gmx_ddbox_t;
struct gmx_mtop_t;
struct t_commrec;
struct t_inputrec;

namespace gmx
{
struct DomdecOptions;
class MDLogger;
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

/*! \brief Checks that requests for PP and PME ranks honor basic expectations
 *
 * Issues a fatal error if there are more PME ranks than PP, or if the
 * count of PP ranks has a prime factor that is too large to be likely
 * to have good performance. */
void checkForValidRankCountRequests(int  numRanksRequested,
                                    bool usingPme,
                                    int  numPmeRanksRequested,
                                    bool checkForLargePrimeFactors);

/*! \brief Return the minimum cell size (in nm) required for DD */
real getDDGridSetupCellSizeLimit(const gmx::MDLogger& mdlog,
                                 bool                 request1DAnd1Pulse,
                                 bool                 bDynLoadBal,
                                 real                 dlb_scale,
                                 const t_inputrec&    ir,
                                 real                 systemInfoCellSizeLimit);

/*! \brief Determines the DD grid setup
 *
 * Either implements settings required by the user, or otherwise
 * chooses estimated optimal number of separate PME ranks and DD grid
 * cell setup, DD cell size limits, and the initial ddbox.
 */
DDGridSetup getDDGridSetup(const gmx::MDLogger&           mdlog,
                           const t_commrec*               cr,
                           int                            numRanksRequested,
                           const gmx::DomdecOptions&      options,
                           const DDSettings&              ddSettings,
                           const DDSystemInfo&            systemInfo,
                           real                           cellSizeLimit,
                           const gmx_mtop_t&              mtop,
                           const t_inputrec&              ir,
                           const matrix                   box,
                           gmx::ArrayRef<const gmx::RVec> xGlobal,
                           gmx_ddbox_t*                   ddbox);

#endif
