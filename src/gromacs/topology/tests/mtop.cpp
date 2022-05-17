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

#include <gtest/gtest.h>

#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"

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

/*! \brief Adds intermolecular bonds, assuming atoms 0 to 5 exist */
void addIntermolecularInteractionBonds(gmx_mtop_t* mtop)
{
    mtop->bIntermolecularInteractions = true;
    mtop->intermolecular_ilist        = std::make_unique<InteractionLists>();
    std::vector<int>& iatoms          = (*mtop->intermolecular_ilist)[F_BONDS].iatoms;
    const int         bondType        = 0;
    iatoms.push_back(bondType);
    iatoms.push_back(0);
    iatoms.push_back(1);
    iatoms.push_back(bondType);
    iatoms.push_back(2);
    iatoms.push_back(3);
    iatoms.push_back(bondType);
    iatoms.push_back(4);
    iatoms.push_back(5);
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
    for (gmx::index i = 0; i < gmx::ssize(atomRanges); i++)
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

TEST(MtopTest, AtomHasPerturbedChargeIn14Interaction)
{
    gmx_moltype_t molt;
    molt.atoms.nr   = 0;
    molt.atoms.atom = nullptr;
    EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(0, molt)) << "empty moltype has none";

    {
        SCOPED_TRACE("Use a moltype that has no perturbed charges at all");
        const int numAtoms = 4;
        init_t_atoms(&molt.atoms, numAtoms, false);
        for (int i = 0; i != numAtoms; ++i)
        {
            molt.atoms.atom[i].q  = 1.0;
            molt.atoms.atom[i].qB = 1.0;
        }
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(0, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(1, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(2, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(3, molt));
    }

    {
        SCOPED_TRACE("Use a moltype that has a perturbed charge, but no interactions");
        molt.atoms.atom[1].q  = 1.0;
        molt.atoms.atom[1].qB = -1.0;
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(0, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(1, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(2, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(3, molt));
    }

    {
        SCOPED_TRACE("Use a moltype that has a perturbed charge but no 1-4 interactions");
        std::array<int, 2> bondAtoms = { 0, 1 };
        molt.ilist[F_BONDS].push_back(0, bondAtoms);
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(0, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(1, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(2, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(3, molt));
    }

    {
        SCOPED_TRACE(
                "Use a moltype that has a perturbed charge, but not on an atom with a 1-4 "
                "interaction");
        std::array<int, 2> pairAtoms = { 2, 3 };
        molt.ilist[F_LJ14].push_back(0, pairAtoms);
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(0, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(1, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(2, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(3, molt));
    }

    {
        SCOPED_TRACE("Use a moltype that has a perturbed charge on an atom with a 1-4 interaction");
        std::array<int, 2> perturbedPairAtoms = { 1, 2 };
        molt.ilist[F_LJ14].push_back(0, perturbedPairAtoms);
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(0, molt));
        EXPECT_TRUE(atomHasPerturbedChargeIn14Interaction(1, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(2, molt));
        EXPECT_FALSE(atomHasPerturbedChargeIn14Interaction(3, molt));
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
        EXPECT_EQ(ilistP.list()[F_BONDS].size(), 0);
        EXPECT_EQ(ilistP.list()[F_SETTLE].size(), 1 * (NRAL(F_SETTLE) + 1));
        count++;
    }
    EXPECT_EQ(count, 1);

    EXPECT_EQ(gmx_mtop_ftype_count(mtop, F_BONDS), 0);
    EXPECT_EQ(gmx_mtop_ftype_count(mtop, F_SETTLE), 2);
    EXPECT_EQ(gmx_mtop_interaction_count(mtop, IF_BOND), 0);
    EXPECT_EQ(gmx_mtop_interaction_count(mtop, IF_CONSTRAINT), 2);
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
            EXPECT_EQ(ilistP.list()[F_BONDS].size(), 0);
            EXPECT_EQ(ilistP.list()[F_SETTLE].size(), 1 * (NRAL(F_SETTLE) + 1));
        }
        else
        {
            EXPECT_EQ(ilistP.nmol(), 1);
            EXPECT_EQ(ilistP.list()[F_BONDS].size(), 3 * (NRAL(F_BONDS) + 1));
            EXPECT_EQ(ilistP.list()[F_SETTLE].size(), 0);
        }
        count++;
    }
    EXPECT_EQ(count, 2);

    EXPECT_EQ(gmx_mtop_ftype_count(mtop, F_BONDS), 3);
    EXPECT_EQ(gmx_mtop_ftype_count(mtop, F_SETTLE), 2);
    EXPECT_EQ(gmx_mtop_interaction_count(mtop, IF_BOND), 3);
    EXPECT_EQ(gmx_mtop_interaction_count(mtop, IF_CONSTRAINT), 2);
    EXPECT_EQ(gmx_mtop_interaction_count(mtop, IF_BOND | IF_CONSTRAINT), 0);
    done_atom(&mtop.moltype[0].atoms);
}

} // namespace

} // namespace test

} // namespace gmx
