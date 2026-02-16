/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * \brief Implementations of LINCS constraint reordering
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "lincs_constraint_group_sizes.h"

#include <cassert>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"


namespace gmx
{

void findConstraintGroupSizes(const int                numConstraints,
                              ArrayRef<const AtomPair> constraintsHost,
                              ArrayRef<int>            constraintGroupSize)
{
    GMX_RELEASE_ASSERT(constraintsHost.size() == constraintGroupSize.size(),
                       "Invalid setup for constraints.");
    const int numConstraintsThreads = constraintsHost.size();
    GMX_RELEASE_ASSERT(numConstraints <= numConstraintsThreads, "Invalid constraint setup");
    /* In order to reduce contention on atomicAdds when running on
     * later AMD devices, we are reordering the constraint array
     * into multiple chunks of constraintGroupSize in order
     * to search adjacent pairs with the same i value.
     *
     * One good example is when you have multiple hydrogen
     * groups (ethane depicted with 2 CH3 hydrogen groups):
     *
     *       H  H
     *       |  |
     *    H--C--C--H
     *       |  |
     *       H  H
     *
     * Carbons are bonded to multiple hydrogen, so each
     * C-H constrain will contend on the C atom 3 times.
     * Atomics performance may wildy vary from device to
     * device, so it's better to reduce reliance on atomics
     * as much as we can.
     *
     * The logic here is to sweep over the constraints
     * and look up hydrogen groups (CH3 groups depicted)
     * and accumulate all updates in the heavy atom before
     * we do an atomicAdd().
     *
     * atomicAdd might still be needed since it's possible
     * that constraints from different thread blocks
     * to update the same heavy atom - need to validate this!
     */
    int      c1           = 1;
    AtomPair previousPair = constraintsHost[0];
    while (c1 < numConstraints)
    {
        AtomPair currentPair       = constraintsHost[c1];
        int      hydrogenGroupSize = 0;
        while ((currentPair.i == previousPair.i) && (currentPair.i != -1))
        {
            hydrogenGroupSize++;
            // if we know we are past the number of constraints, break out of the
            // loop here and continue below. Catches edge case of constraints
            // fitting exactly into the buffer
            if ((c1 + hydrogenGroupSize) >= numConstraintsThreads)
            {
                break;
            }
            currentPair = constraintsHost[c1 + hydrogenGroupSize];
        }
        previousPair = currentPair;
        // now advances c1 to the first atom in the next hydrogen group
        // and set the size of the group in the array. All other entries
        // are set to zero.
        // If no constraint group is present, then the sentinel value of
        // -1 stays unchanged in the array to indicate that no group has
        // been found.
        if (hydrogenGroupSize > 0)
        {
            // Index into c1 - 1, as the very first constraint by definition
            // can't match any previous constraints
            constraintGroupSize[c1 - 1] = hydrogenGroupSize;
            GMX_RELEASE_ASSERT(c1 + hydrogenGroupSize - 1 < numConstraintsThreads,
                               "Size of constraint group overflows storage");
            for (int cgc = 0; cgc < hydrogenGroupSize; cgc++)
            {
                constraintGroupSize[c1 + cgc] = 0;
            }
            c1 += hydrogenGroupSize;
        }
        c1++;
    }
}

} // namespace gmx
