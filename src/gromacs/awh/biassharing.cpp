/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief
 * Implements bias sharing checking functionality.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "biassharing.h"

#include <vector>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{

bool haveBiasSharingWithinSimulation(const AwhParams &awhParams)
{
    bool haveSharing = false;

    for (int k = 0; k < awhParams.numBias; k++)
    { 
        int shareGroup = awhParams.awhBiasParams[k].shareGroup;
        if (shareGroup > 0)
        {
            for (int i = k + 1; i < awhParams.numBias; i++)
            {
                if (awhParams.awhBiasParams[i].shareGroup == shareGroup)
                {
                    haveSharing = true;
                }
            }
        }
    }

    return haveSharing;
}

void checkBiasSharingMultiSim(const AwhParams           &awhParams,
                              const std::vector<size_t> &pointSize,
                              const gmx_multisim_t      *multiSimComm)
{
    const int        numSim = multiSimComm->nsim;

    std::vector<int> share(awhParams.numBias);
    for (int b = 0; b < awhParams.numBias; b++)
    {
        share[b] = (awhParams.awhBiasParams[b].shareGroup > 0 ? 1 : 0);
    }
    gmx_sumi_sim(share.size(), share.data(), multiSimComm);
    for (int b = 0; b < awhParams.numBias; b++)
    {
        if (share[b] > 0 && share[b] != multiSimComm->nsim)
        {
            gmx_fatal(FARGS, "Some, but not all, of the %d simulations share AWH bias %d. Either no or all simulations should share a bias.\n", numSim, b + 1);
        }
    }

    std::vector<int> intervals(numSim*2);
    intervals[numSim*0 + multiSimComm->sim] = awhParams.nstSampleCoord;
    intervals[numSim*1 + multiSimComm->sim] = awhParams.numSamplesUpdateFreeEnergy;
    gmx_sumi_sim(intervals.size(), intervals.data(), multiSimComm);
    for (int sim = 1; sim < numSim; sim++)
    {
        if (intervals[sim] != intervals[0])
        {
            gmx_fatal(FARGS, "All simulations should have the same AWH sample interval");
        }
        if (intervals[numSim + sim] != intervals[numSim])
        {
            gmx_fatal(FARGS, "All simulations should have the same AWH free-energy update interval");
        }
    }

    /* Check the point sizes. This is a sufficient condition for running
     * as shared multi-sim run. No physics checks are performed here.
     */
    for (int b = 0; b < awhParams.numBias; b++)
    {
        if (awhParams.awhBiasParams[b].shareGroup > 0)
        {
            std::vector<gmx_int64_t> pointSizes(numSim);
            pointSizes[multiSimComm->sim] = pointSize[b];
            gmx_sumli_sim(pointSizes.size(), pointSizes.data(), multiSimComm);
            for (int sim = 1; sim < numSim; sim++)
            {
                if (pointSizes[sim] != pointSizes[0])
                {
                    gmx_fatal(FARGS, "Shared AWH bias %d has different grid sizes in different simulations\n", b + 1);
                }
            }
        }
    }
}

} // namespace gmx
