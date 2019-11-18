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
 * \brief This file declares command-line options for mdrun related to
 * domain decomposition.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_OPTIONS_H
#define GMX_DOMDEC_OPTIONS_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \brief The options for the domain decomposition MPI task ordering. */
enum class DdRankOrder
{
    select,     //!< First value (needed to cope with command-line parsing)
    interleave, //!< Interleave the PP and PME ranks
    pp_pme,     //!< First all PP ranks, all PME rank at the end
    cartesian,  //!< Use Cartesian communicators for PP, PME and PP-PME
    Count       //!< The number of options
};

/*! \brief The options for the dynamic load balancing. */
enum class DlbOption
{
    select,           //!< First value (needed to cope with command-line parsing)
    turnOnWhenUseful, //!< Turn on DLB when we think it would improve performance
    no,               //!< Never turn on DLB
    yes,              //!< Turn on DLB from the start and keep it on
    Count             //!< The number of options
};

/*! \libinternal \brief Structure containing all (command line) options for the domain decomposition */
struct DomdecOptions
{
    //! If true, check that all bonded interactions have been assigned to exactly one domain/rank.
    bool checkBondedInteractions = true;
    //! If true, don't communicate all atoms between the non-bonded cut-off and the larger bonded cut-off, but only those that have non-local bonded interactions. This significantly reduces the communication volume.
    bool useBondedCommunication = true;
    //! The domain decomposition grid cell count, 0 means let domdec choose based on the number of ranks.
    ivec numCells = { 0 };
    //! The number of separate PME ranks requested, -1 = auto.
    int numPmeRanks = -1;
    //! Ordering of the PP and PME ranks, values from enum above.
    DdRankOrder rankOrder = DdRankOrder::interleave;
    //! The minimum communication range, used for extended the communication range for bonded interactions (nm).
    real minimumCommunicationRange = 0;
    //! Communication range for atom involved in constraints (P-LINCS) (nm).
    real constraintCommunicationRange = 0;
    //! Dynamic load balancing option, values from enum above.
    DlbOption dlbOption = DlbOption::turnOnWhenUseful;
    /*! \brief Fraction in (0,1) by whose reciprocal the initial
     * DD cell size will be increased in order to provide a margin
     * in which dynamic load balancing can act, while preserving
     * the minimum cell size. */
    real dlbScaling = 0.8;
    //! String containing a vector of the relative sizes in the x direction of the corresponding DD cells.
    const char* cellSizeX = nullptr;
    //! String containing a vector of the relative sizes in the y direction of the corresponding DD cells.
    const char* cellSizeY = nullptr;
    //! String containing a vector of the relative sizes in the z direction of the corresponding DD cells.
    const char* cellSizeZ = nullptr;
};

} // namespace gmx

#endif
