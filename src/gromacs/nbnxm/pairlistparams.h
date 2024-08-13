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
 * \brief
 * Declares the PairlistType enum and PairlistParams class
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_PAIRLISTPARAMS_H
#define GMX_NBNXM_PAIRLISTPARAMS_H

#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

#include "nbnxm_enums.h"

namespace gmx
{

enum class NbnxmKernelType;

//! Gives the i-cluster size for each pairlist type
static constexpr gmx::EnumerationArray<PairlistType, int> IClusterSizePerListType = {
    { 4, 4, 4, c_nbnxnGpuClusterSize }
};
//! Gives the j-cluster size for each pairlist type
static constexpr gmx::EnumerationArray<PairlistType, int> JClusterSizePerListType = {
    { 2, 4, 8, c_nbnxnGpuClusterSize }
};
//! True if given pairlist type is used on GPU, false if on CPU.
static constexpr gmx::EnumerationArray<PairlistType, bool> sc_isGpuPairListType = {
    { false, false, false, true }
};

/*! \internal
 * \brief The setup for generating and pruning the nbnxn pair list.
 *
 * Without dynamic pruning rlistOuter=rlistInner.
 */
struct PairlistParams
{
    /*! \brief Constructor producing a struct with dynamic pruning disabled
     */
    PairlistParams(NbnxmKernelType kernelType, bool haveFep, real rlist, bool haveMultipleDomains);

    //! The type of cluster-pair list
    PairlistType pairlistType;
    //! Tells whether we have perturbed interactions
    bool haveFep_;
    //! Cut-off of the larger, outer pair-list
    real rlistOuter;
    //! Cut-off of the smaller, inner pair-list
    real rlistInner;
    //! True when using DD with multiple domains
    bool haveMultipleDomains_;
    //! Are we using dynamic pair-list pruning
    bool useDynamicPruning;
    //! The interval in steps for computing non-bonded interactions, =1 without MTS
    int mtsFactor;
    //! Pair-list dynamic pruning interval
    int nstlistPrune;
    //! The number parts to divide the pair-list into for rolling pruning, a value of 1 gives no rolling pruning
    int numRollingPruningParts;
    //! Lifetime in steps of the pair-list
    int lifetime;
};

} // namespace gmx

#endif
