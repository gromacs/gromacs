/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2008- The GROMACS Authors
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
#include "gmxpre.h"

#include "constraint_gpu_helpers.h"

#include <array>
#include <filesystem>
#include <utility>

#include "gromacs/topology/atoms.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"

int countCoupled(int                                                a,
                 gmx::ArrayRef<int>                                 numCoupledConstraints,
                 const gmx::ListOfLists<AtomsAdjacencyListElement>& atomsAdjacencyList)

{
    int counted = 0;
    for (const auto& adjacentAtom : atomsAdjacencyList[a])
    {
        const int c2 = adjacentAtom.indexOfConstraint_;
        if (numCoupledConstraints[c2] == -1)
        {
            numCoupledConstraints[c2] = 0; // To indicate we've been here
            counted += 1
                       + countCoupled(adjacentAtom.indexOfSecondConstrainedAtom_,
                                      numCoupledConstraints,
                                      atomsAdjacencyList);
        }
    }
    return counted;
}

std::vector<int> countNumCoupledConstraints(gmx::ArrayRef<const int> iatoms,
                                            const gmx::ListOfLists<AtomsAdjacencyListElement>& atomsAdjacencyList)
{
    const int        stride         = 1 + NRAL(F_CONSTR);
    const int        numConstraints = iatoms.ssize() / stride;
    std::vector<int> numCoupledConstraints(numConstraints, -1);
    for (int c = 0; c < numConstraints; c++)
    {
        const int a1 = iatoms[stride * c + 1];
        const int a2 = iatoms[stride * c + 2];
        if (numCoupledConstraints[c] == -1)
        {
            numCoupledConstraints[c] = countCoupled(a1, numCoupledConstraints, atomsAdjacencyList)
                                       + countCoupled(a2, numCoupledConstraints, atomsAdjacencyList);
        }
    }

    return numCoupledConstraints;
}

//! Constructs and returns an atom constraint adjacency list
gmx::ListOfLists<AtomsAdjacencyListElement> constructAtomsAdjacencyList(const int numAtoms,
                                                                        gmx::ArrayRef<const int> iatoms)
{
    const int stride         = 1 + NRAL(F_CONSTR);
    const int numConstraints = iatoms.ssize() / stride;

    // count how many constraints each atom has. These counts will be used to place the constraints
    // for each atom in contiguous memory in a ListOfLists object, which is more performant than
    // using vector<vector<>>
    std::vector<int> count(numAtoms + 1);
    count[0] = 0;
    for (int c = 0; c < numConstraints; c++)
    {
        int a1 = iatoms[stride * c + 1];
        int a2 = iatoms[stride * c + 2];
        count[a1 + 1]++;
        count[a2 + 1]++;
    }

    // perform an inclusive prefix sum on count to use for tracking insertion indices
    std::inclusive_scan(count.begin(), count.end(), count.begin());
    // take a copy of the initial state of count as these are the list ranges
    std::vector<int> listRanges = count;
    // initialise the elements list. The dummy value of (0,0,0) is required as AtomsAdjacencyListElement
    // has no default constructor, but all these elements should be overwritten in the loop below
    std::vector<AtomsAdjacencyListElement> elements(listRanges.back(), AtomsAdjacencyListElement(0, 0, 0));

    for (int c = 0; c < numConstraints; c++)
    {
        int a1 = iatoms[stride * c + 1];
        int a2 = iatoms[stride * c + 2];

        // Each constraint will be represented as a tuple, containing index of the second
        // constrained atom, index of the constraint and a sign that indicates the order of atoms in
        // which they are listed. Sign is needed to compute the mass factors.
        elements[count[a1]++] = AtomsAdjacencyListElement(a2, c, +1);
        elements[count[a2]++] = AtomsAdjacencyListElement(a1, c, -1);
    }
    return gmx::ListOfLists<AtomsAdjacencyListElement>(std::move(listRanges), std::move(elements));
}

bool isNumCoupledConstraintsSupported(const gmx_mtop_t& mtop, int threadsPerBlock)
{
    for (const gmx_moltype_t& molType : mtop.moltype)
    {
        gmx::ArrayRef<const int> iatoms = molType.ilist[F_CONSTR].iatoms;
        const auto atomsAdjacencyList   = constructAtomsAdjacencyList(molType.atoms.nr, iatoms);
        // Compute, how many constraints are coupled to each constraint
        const auto numCoupledConstraints = countNumCoupledConstraints(iatoms, atomsAdjacencyList);
        for (const int numCoupled : numCoupledConstraints)
        {
            if (numCoupled > threadsPerBlock)
            {
                return false;
            }
        }
    }

    return true;
}

int computeTotalNumSettles(const gmx_mtop_t& mtop)
{ // This is to prevent the assertion failure for the systems without water
    int totalSettles = 0;
    for (unsigned mt = 0; mt < mtop.moltype.size(); mt++)
    {
        const int        nral1           = 1 + NRAL(F_SETTLE);
        InteractionList  interactionList = mtop.moltype[mt].ilist[F_SETTLE];
        std::vector<int> iatoms          = interactionList.iatoms;
        totalSettles += iatoms.size() / nral1;
    }
    return totalSettles;
}

SettleWaterTopology getSettleTopologyData(const gmx_mtop_t& mtop)
{
    real mO = -1.0;
    real mH = -1.0;

    for (unsigned mt = 0; mt < mtop.moltype.size(); mt++)
    {
        const int        nral1           = 1 + NRAL(F_SETTLE);
        InteractionList  interactionList = mtop.moltype[mt].ilist[F_SETTLE];
        std::vector<int> iatoms          = interactionList.iatoms;
        for (unsigned i = 0; i < iatoms.size() / nral1; i++)
        {
            WaterMolecule settler;
            settler.ow1 = iatoms[i * nral1 + 1]; // Oxygen index
            settler.hw2 = iatoms[i * nral1 + 2]; // First hydrogen index
            settler.hw3 = iatoms[i * nral1 + 3]; // Second hydrogen index
            t_atom ow1  = mtop.moltype[mt].atoms.atom[settler.ow1];
            t_atom hw2  = mtop.moltype[mt].atoms.atom[settler.hw2];
            t_atom hw3  = mtop.moltype[mt].atoms.atom[settler.hw3];

            if (mO < 0)
            {
                mO = ow1.m;
            }
            if (mH < 0)
            {
                mH = hw2.m;
            }
            GMX_RELEASE_ASSERT(mO == ow1.m,
                               "Topology has different values for oxygen mass. Should be identical "
                               "in order to use SETTLE.");
            GMX_RELEASE_ASSERT(hw2.m == hw3.m && hw2.m == mH,
                               "Topology has different values for hydrogen mass. Should be "
                               "identical in order to use SETTLE.");
        }
    }
    GMX_RELEASE_ASSERT(mO > 0, "Could not find oxygen mass in the topology. Needed in SETTLE.");
    GMX_RELEASE_ASSERT(mH > 0, "Could not find hydrogen mass in the topology. Needed in SETTLE.");

    // TODO Very similar to SETTLE initialization on CPU. Should be lifted to a separate method
    // (one that gets dOH and dHH values and checks them for consistency)
    int settle_type = -1;
    for (unsigned mt = 0; mt < mtop.moltype.size(); mt++)
    {
        const int       nral1           = 1 + NRAL(F_SETTLE);
        InteractionList interactionList = mtop.moltype[mt].ilist[F_SETTLE];
        for (int i = 0; i < interactionList.size(); i += nral1)
        {
            if (settle_type == -1)
            {
                settle_type = interactionList.iatoms[i];
            }
            else if (interactionList.iatoms[i] != settle_type)
            {
                gmx_fatal(FARGS,
                          "The [molecules] section of your topology specifies more than one block "
                          "of\n"
                          "a [moleculetype] with a [settles] block. Only one such is allowed.\n"
                          "If you are trying to partition your solvent into different *groups*\n"
                          "(e.g. for freezing, T-coupling, etc.), you are using the wrong "
                          "approach. Index\n"
                          "files specify groups. Otherwise, you may wish to change the least-used\n"
                          "block of molecules with SETTLE constraints into 3 normal constraints.");
            }
        }
    }

    GMX_RELEASE_ASSERT(settle_type >= 0, "settle_init called without settles");

    real dOH = mtop.ffparams.iparams[settle_type].settle.doh;
    real dHH = mtop.ffparams.iparams[settle_type].settle.dhh;

    return { mO, mH, dOH, dHH };
}

LocalSettleData computeNumSettles(const InteractionDefinitions& idef)
{
    const int              nral1     = 1 + NRAL(F_SETTLE);
    const InteractionList& il_settle = idef.il[F_SETTLE];
    return { il_settle.size() / nral1, nral1 };
}

gmx::ArrayRef<const int> localSettleAtoms(const InteractionDefinitions& idef)
{
    const InteractionList& il_settle = idef.il[F_SETTLE];
    return il_settle.iatoms;
}
