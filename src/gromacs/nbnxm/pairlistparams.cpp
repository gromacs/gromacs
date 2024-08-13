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
 * \brief
 * Implements the PairlistParams constructor
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#include "gmxpre.h"

#include "pairlistparams.h"

#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/utility/gmxassert.h"

#include "nbnxm_geometry.h"

namespace gmx
{

PairlistParams::PairlistParams(const NbnxmKernelType kernelType,
                               const bool            haveFep,
                               const real            rlist,
                               const bool            haveMultipleDomains) :
    haveFep_(haveFep),
    rlistOuter(rlist),
    rlistInner(rlist),
    haveMultipleDomains_(haveMultipleDomains),
    useDynamicPruning(false),
    mtsFactor(1),
    nstlistPrune(-1),
    numRollingPruningParts(1),
    lifetime(-1)
{
    if (!kernelTypeUsesSimplePairlist(kernelType))
    {
        pairlistType = PairlistType::HierarchicalNxN;
    }
    else
    {
        switch (sc_jClusterSize(kernelType))
        {
            case 2: pairlistType = PairlistType::Simple4x2; break;
            case 4: pairlistType = PairlistType::Simple4x4; break;
            case 8: pairlistType = PairlistType::Simple4x8; break;
            default: GMX_RELEASE_ASSERT(false, "Kernel type does not have a pairlist type");
        }
    }
}

} // namespace gmx
