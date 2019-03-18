/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief
 * Declares the PairlistParams class
 *
 * This class holds the Nbnxm pairlist parameters.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_PAIRLISTPARAMS_H
#define GMX_NBNXM_PAIRLISTPARAMS_H

#include "gromacs/utility/real.h"

#include "locality.h"

enum class PairlistType;

namespace Nbnxm
{
enum class KernelType;
}


/*! \internal
 * \brief The setup for generating and pruning the nbnxn pair list.
 *
 * Without dynamic pruning rlistOuter=rlistInner.
 */
struct PairlistParams
{
    /*! \brief Constructor producing a struct with dynamic pruning disabled
     */
    PairlistParams(Nbnxm::KernelType kernelType,
                   bool              haveFep,
                   real              rlist,
                   bool              haveMultipleDomains);

    PairlistType pairlistType;           //!< The type of cluster-pair list
    bool         haveFep;                //!< Tells whether we have perturbed interactions
    real         rlistOuter;             //!< Cut-off of the larger, outer pair-list
    real         rlistInner;             //!< Cut-off of the smaller, inner pair-list
    bool         haveMultipleDomains;    //!< True when using DD with multiple domains
    bool         useDynamicPruning;      //!< Are we using dynamic pair-list pruning
    int          nstlistPrune;           //!< Pair-list dynamic pruning interval
    int          numRollingPruningParts; //!< The number parts to divide the pair-list into for rolling pruning, a value of 1 gives no rolling pruning
    int          lifetime;               //!< Lifetime in steps of the pair-list
};

#endif
