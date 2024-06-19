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

/*! \brief Add an ethane molecule with uniform charges on all atoms and states */
void addEthane(gmx_mtop_t* mtop)
{
    t_iparams iparams;
    // Add a constraint type
    iparams.constr = { 0.5, 0.5 };
    mtop->ffparams.iparams.push_back(iparams);
    const int constraintInteractionType = 0;
    // Add an angle type
    iparams.harmonic = { 120, 1000, 120, 1000 };
    mtop->ffparams.iparams.push_back(iparams);
    const int angleInteractionType = 1;
    // Add an LJ-14 type
    iparams.lj14 = { 1, 100, 1, 100 };
    mtop->ffparams.iparams.push_back(iparams);
    const int oneFourInteractionType = 2;

    const int numAtoms = 8;
    const int carbon0  = 0;
    const int carbon1  = 4;
    // Atoms 1,2,3,5,6,7 are hydrogens.

    gmx_moltype_t& ethaneMoleculeType = mtop->moltype.emplace_back(gmx_moltype_t{});

    // Describe the interactions
    // clang-format off
    ethaneMoleculeType.ilist[F_CONSTR].iatoms = { constraintInteractionType, carbon0, 1,
                                                  constraintInteractionType, carbon0, 2,
                                                  constraintInteractionType, carbon0, 3,
                                                  constraintInteractionType, carbon1, 5,
                                                  constraintInteractionType, carbon1, 6,
                                                  constraintInteractionType, carbon1, 7 };
    ethaneMoleculeType.ilist[F_ANGLES].iatoms = { angleInteractionType, 1, carbon0, 2,
                                                  angleInteractionType, 1, carbon0, 3,
                                                  angleInteractionType, 2, carbon0, 3,
                                                  angleInteractionType, 5, carbon1, 6,
                                                  angleInteractionType, 5, carbon1, 7,
                                                  angleInteractionType, 6, carbon1, 7 };
    ethaneMoleculeType.ilist[F_LJ14].iatoms = { oneFourInteractionType, 1, 5,
                                                oneFourInteractionType, 1, 6,
                                                oneFourInteractionType, 1, 7,
                                                oneFourInteractionType, 2, 5,
                                                oneFourInteractionType, 2, 6,
                                                oneFourInteractionType, 2, 7,
                                                oneFourInteractionType, 3, 5,
                                                oneFourInteractionType, 3, 6,
                                                oneFourInteractionType, 3, 7 };
    // clang-format on

    // Set up the atoms, including arbitrary uniform charges
    ethaneMoleculeType.atoms.nr = numAtoms;
    snew(ethaneMoleculeType.atoms.atom, numAtoms);
    for (int i = 0; i != numAtoms; ++i)
    {
        ethaneMoleculeType.atoms.atom[i].q  = 1.0;
        ethaneMoleculeType.atoms.atom[i].qB = 1.0;
    }
}

TEST(MtopTest, CanSortPerturbedInteractionsCorrectly)
{
    gmx_mtop_t mtop;

    // Add a molecule type of normal ethane
    addEthane(&mtop);
    // Perturb the charge on atoms 1 and 7 by changing the B-state
    // value.
    const auto& ethaneMoleculeType = mtop.moltype.back();
    ethaneMoleculeType.atoms.atom[1].qB *= -1.0;
    ethaneMoleculeType.atoms.atom[7].qB *= -1.0;

    // Add a molecule block of ethane with perturbed charge
    gmx_molblock_t molblock;
    molblock.type = 0;
    molblock.nmol = 1;
    mtop.molblock.push_back(molblock);
    mtop.natoms += mtop.moltype[molblock.type].atoms.nr;
    mtop.finalize();

    SCOPED_TRACE("Make a local topology without sorting perturbed interactions to the end");
    {
        gmx_localtop_t localTopology(mtop.ffparams);
        const bool     perturbedInteractionsAtEnd = false;
        gmx_mtop_generate_local_top(mtop, &localTopology, perturbedInteractionsAtEnd);
        ASSERT_EQ(27, localTopology.idef.il[F_LJ14].size());
        int iatomsIndex = 0;
        // All the interactions remain in the order in which they were
        // inserted.
        for (const int i : { 1, 2, 3 })
        {
            for (const int j : { 5, 6, 7 })
            {
                ++iatomsIndex;
                EXPECT_EQ(i, localTopology.idef.il[F_LJ14].iatoms[iatomsIndex++]);
                EXPECT_EQ(j, localTopology.idef.il[F_LJ14].iatoms[iatomsIndex++]);
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
        ASSERT_EQ(27, localTopology.idef.il[F_LJ14].size());
        int iatomsIndex = 0;
        // All the interactions including neither atom 1 or 7 remain in
        // the order in which they were inserted.
        for (const int i : { 2, 3 })
        {
            for (const int j : { 5, 6 })
            {
                ++iatomsIndex;
                EXPECT_EQ(i, localTopology.idef.il[F_LJ14].iatoms[iatomsIndex++]);
                EXPECT_EQ(j, localTopology.idef.il[F_LJ14].iatoms[iatomsIndex++]);
            }
        }
        // Then all interactions with atom 1
        for (const int i : { 1 })
        {
            for (const int j : { 5, 6, 7 })
            {
                ++iatomsIndex;
                EXPECT_EQ(i, localTopology.idef.il[F_LJ14].iatoms[iatomsIndex++]);
                EXPECT_EQ(j, localTopology.idef.il[F_LJ14].iatoms[iatomsIndex++]);
            }
        }
        // Then all interactions with atom 7
        for (const int i : { 2, 3 })
        {
            for (const int j : { 7 })
            {
                ++iatomsIndex;
                EXPECT_EQ(i, localTopology.idef.il[F_LJ14].iatoms[iatomsIndex++]);
                EXPECT_EQ(j, localTopology.idef.il[F_LJ14].iatoms[iatomsIndex++]);
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
