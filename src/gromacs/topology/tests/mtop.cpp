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
 * Implements test of mtop routines
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_topology
 */
#include "gmxpre.h"

#include <array>
#include <initializer_list>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/topology/atoms.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/topologyhelpers.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief
 * Creates dummy topology with two differently sized residues.
 *
 * Residue begin and end are set to allow checking routines
 * that make use of the boundaries.
 *
 * \returns The residue ranges.
 */
std::vector<gmx::Range<int>> createTwoResidueTopology(gmx_mtop_t* mtop)
{
    auto& moltype        = mtop->moltype.emplace_back();
    int   residueOneSize = 5;
    int   residueTwoSize = 4;
    moltype.atoms.nr     = residueOneSize + residueTwoSize;
    snew(moltype.atoms.atom, residueOneSize + residueTwoSize);
    for (int i = 0; i < residueOneSize; i++)
    {
        moltype.atoms.atom[i].resind = 0;
    }
    for (int i = residueOneSize; i < residueOneSize + residueTwoSize; i++)
    {
        moltype.atoms.atom[i].resind = 1;
    }

    mtop->molblock.resize(1);
    mtop->molblock[0].type = 0;
    mtop->molblock[0].nmol = 1;
    mtop->natoms           = moltype.atoms.nr * mtop->molblock[0].nmol;
    std::vector<gmx::Range<int>> residueRange;
    residueRange.emplace_back(0, residueOneSize);
    residueRange.emplace_back(residueOneSize, residueOneSize + residueTwoSize);
    return residueRange;
}

TEST(MtopTest, RangeBasedLoop)
{
    gmx_mtop_t mtop;
    addNWaterMolecules(&mtop, 2);
    mtop.finalize();
    int count = 0;
    for (const AtomProxy atomP : AtomRange(mtop))
    {
        EXPECT_EQ(atomP.globalAtomNumber(), count);
        ++count;
    }
    EXPECT_EQ(count, 6);
    done_atom(&mtop.moltype[0].atoms);
}

TEST(MtopTest, Operators)
{
    gmx_mtop_t mtop;
    addNWaterMolecules(&mtop, 2);
    mtop.finalize();
    AtomIterator it(mtop);
    AtomIterator otherIt(mtop);
    EXPECT_EQ((*it).globalAtomNumber(), 0);
    EXPECT_EQ(it->globalAtomNumber(), 0);
    EXPECT_TRUE(it == otherIt);
    EXPECT_FALSE(it != otherIt);
    ++it;
    EXPECT_EQ(it->globalAtomNumber(), 1);
    it++;
    EXPECT_EQ(it->globalAtomNumber(), 2);
    EXPECT_TRUE(it != otherIt);
    EXPECT_FALSE(it == otherIt);
    done_atom(&mtop.moltype[0].atoms);
}

TEST(MtopTest, CanFindResidueStartAndEndAtoms)
{
    gmx_mtop_t mtop;
    auto       expectedResidueRange = createTwoResidueTopology(&mtop);
    mtop.finalize();

    auto atomRanges = atomRangeOfEachResidue(mtop.moltype[0]);
    ASSERT_EQ(atomRanges.size(), expectedResidueRange.size());
    for (gmx::Index i = 0; i < gmx::ssize(atomRanges); i++)
    {
        EXPECT_EQ(atomRanges[i].begin(), expectedResidueRange[i].begin());
        EXPECT_EQ(atomRanges[i].end(), expectedResidueRange[i].end());
        ASSERT_EQ(atomRanges[i].size(), expectedResidueRange[i].size());
    }
    int rangeIndex = 0;
    for (const auto& range : atomRanges)
    {
        auto referenceRangeIt = expectedResidueRange[rangeIndex].begin();
        for (int i : range)
        {
            EXPECT_EQ(i, *referenceRangeIt);
            referenceRangeIt++;
        }
        rangeIndex++;
    }
}

TEST(MtopTest, CanSortPerturbedInteractionsCorrectly)
{
    gmx_mtop_t mtop;

    // Add a molecule type of normal ethane
    const bool bondsAreFlexible = false;
    addEthane(&mtop, bondsAreFlexible);
    // Perturb the charge on atoms 1 and 7 by changing the B-state
    // value.
    const auto& ethaneMoleculeType = mtop.moltype.back();
    ethaneMoleculeType.atoms.atom[1].qB *= -1.0;
    ethaneMoleculeType.atoms.atom[7].qB *= -1.0;
    mtop.finalize();

    SCOPED_TRACE("Make a local topology without sorting perturbed interactions to the end");
    {
        gmx_localtop_t localTopology(mtop.ffparams);
        const bool     perturbedInteractionsAtEnd = false;
        gmx_mtop_generate_local_top(mtop, &localTopology, perturbedInteractionsAtEnd);
        ASSERT_EQ(27, localTopology.idef.il[InteractionFunction::LennardJones14].size());
        int iatomsIndex = 0;
        // All the interactions remain in the order in which they were
        // inserted.
        for (const int i : { 1, 2, 3 })
        {
            for (const int j : { 5, 6, 7 })
            {
                ++iatomsIndex;
                EXPECT_EQ(i, localTopology.idef.il[InteractionFunction::LennardJones14].iatoms[iatomsIndex++]);
                EXPECT_EQ(j, localTopology.idef.il[InteractionFunction::LennardJones14].iatoms[iatomsIndex++]);
            }
        }
    }

    SCOPED_TRACE(
            "Make a local topology that sorts perturbed interactions to the end, note that the "
            "interactions including atoms 1 or 7 are now at the end");
    {
        gmx_localtop_t localTopology(mtop.ffparams);
        const bool     perturbedInteractionsAtEnd = true;
        gmx_mtop_generate_local_top(mtop, &localTopology, perturbedInteractionsAtEnd);
        ASSERT_EQ(27, localTopology.idef.il[InteractionFunction::LennardJones14].size());
        int iatomsIndex = 0;
        // All the interactions including neither atom 1 or 7 remain in
        // the order in which they were inserted.
        for (const int i : { 2, 3 })
        {
            for (const int j : { 5, 6 })
            {
                ++iatomsIndex;
                EXPECT_EQ(i, localTopology.idef.il[InteractionFunction::LennardJones14].iatoms[iatomsIndex++]);
                EXPECT_EQ(j, localTopology.idef.il[InteractionFunction::LennardJones14].iatoms[iatomsIndex++]);
            }
        }
        // Then all interactions with atom 1
        for (const int i : { 1 })
        {
            for (const int j : { 5, 6, 7 })
            {
                ++iatomsIndex;
                EXPECT_EQ(i, localTopology.idef.il[InteractionFunction::LennardJones14].iatoms[iatomsIndex++]);
                EXPECT_EQ(j, localTopology.idef.il[InteractionFunction::LennardJones14].iatoms[iatomsIndex++]);
            }
        }
        // Then all interactions with atom 7
        for (const int i : { 2, 3 })
        {
            for (const int j : { 7 })
            {
                ++iatomsIndex;
                EXPECT_EQ(i, localTopology.idef.il[InteractionFunction::LennardJones14].iatoms[iatomsIndex++]);
                EXPECT_EQ(j, localTopology.idef.il[InteractionFunction::LennardJones14].iatoms[iatomsIndex++]);
            }
        }
    }
}

TEST(IListRangeTest, RangeBasedLoopWorks)
{
    gmx_mtop_t mtop;
    addNWaterMolecules(&mtop, 2);
    mtop.finalize();

    int count = 0;
    for (const IListProxy ilistP : IListRange(mtop))
    {
        EXPECT_EQ(ilistP.nmol(), 2);
        EXPECT_EQ(ilistP.list()[InteractionFunction::Bonds].size(), 0);
        EXPECT_EQ(ilistP.list()[InteractionFunction::SETTLE].size(),
                  1 * (NRAL(InteractionFunction::SETTLE) + 1));
        count++;
    }
    EXPECT_EQ(count, 1);

    EXPECT_EQ(gmx_mtop_ftype_count(mtop, InteractionFunction::Bonds), 0);
    EXPECT_EQ(gmx_mtop_ftype_count(mtop, InteractionFunction::SETTLE), 2);
    EXPECT_EQ(gmx_mtop_interaction_count(mtop, sc_ifcBondedForce), 0);
    EXPECT_EQ(gmx_mtop_interaction_count(mtop, sc_ifcConstraint), 2);
    done_atom(&mtop.moltype[0].atoms);
}

TEST(IListRangeTest, RangeBasedLoopWithIntermolecularInteraction)
{
    gmx_mtop_t mtop;
    addNWaterMolecules(&mtop, 2);
    addIntermolecularInteractionBonds(&mtop);
    mtop.finalize();

    int count = 0;
    for (const IListProxy ilistP : IListRange(mtop))
    {
        if (count == 0)
        {
            EXPECT_EQ(ilistP.nmol(), 2);
            EXPECT_EQ(ilistP.list()[InteractionFunction::Bonds].size(), 0);
            EXPECT_EQ(ilistP.list()[InteractionFunction::SETTLE].size(),
                      1 * (NRAL(InteractionFunction::SETTLE) + 1));
        }
        else
        {
            EXPECT_EQ(ilistP.nmol(), 1);
            EXPECT_EQ(ilistP.list()[InteractionFunction::Bonds].size(),
                      3 * (NRAL(InteractionFunction::Bonds) + 1));
            EXPECT_EQ(ilistP.list()[InteractionFunction::SETTLE].size(), 0);
        }
        count++;
    }
    EXPECT_EQ(count, 2);

    EXPECT_EQ(gmx_mtop_ftype_count(mtop, InteractionFunction::Bonds), 3);
    EXPECT_EQ(gmx_mtop_ftype_count(mtop, InteractionFunction::SETTLE), 2);
    EXPECT_EQ(gmx_mtop_interaction_count(mtop, sc_ifcBondedForce), 3);
    EXPECT_EQ(gmx_mtop_interaction_count(mtop, sc_ifcConstraint), 2);
    EXPECT_EQ(gmx_mtop_interaction_count(mtop, sc_ifcBondedForce | sc_ifcConstraint), 0);
    done_atom(&mtop.moltype[0].atoms);
}

} // namespace

} // namespace test

} // namespace gmx
