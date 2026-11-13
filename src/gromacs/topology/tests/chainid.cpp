/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Tests chain ID construction routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_topology
 */
#include "gmxpre.h"

#include "gromacs/topology/chainid.h"

#include <gtest/gtest.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"

#include "testutils/testfilemanager.h"
#include "testutils/tprfilegenerator.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Test helper compound type for a sequence of particles with the same chain id
struct ParticleChain
{
    //! Length of sequence
    int length;
    //! Chain ID for particles in sequence
    char id;
};

//! \brief Test fixture for functions constructing chain IDs from topology data
class ChainIdFromTopologyTest : public ::testing::Test
{
public:
    //! Test all chain ID construction functions for topology data read from \p inputConf
    static void runTest(const std::filesystem::path& inputConf, ArrayRef<const ParticleChain> expectedChains)
    {
        bool       haveFullTopology;
        gmx_mtop_t mtop;
        t_topology top;
        PbcType    pbcType;
        matrix     box;

        readConfAndTopology(inputConf, &haveFullTopology, &mtop, &pbcType, nullptr, nullptr, box);
        read_tps_conf(inputConf, &top, &pbcType, nullptr, nullptr, box, false);

        runMtopLegacyTest(mtop, expectedChains);
        runMtopModernTest(mtop, expectedChains);
        runTopologyTest(top, expectedChains);

        done_top(&top);
    }

private:
    //! Test makeChainIdentifiersAfterTprReading for input mtop
    static void runMtopModernTest(const gmx_mtop_t& mtop, ArrayRef<const ParticleChain> expectedChains)
    {
        t_atoms           atoms = gmx_mtop_global_atoms(mtop);
        RangePartitioning mols  = gmx_mtop_molecules(mtop);

        makeChainIdentifiersAfterTprReading(&atoms, mols);
        checkChainIds(atoms, expectedChains);

        done_atom(&atoms);
    }

    //! Test tpx_make_chain_identifiers for input mtop
    static void runMtopLegacyTest(const gmx_mtop_t& mtop, ArrayRef<const ParticleChain> expectedChains)
    {
        t_atoms atoms = gmx_mtop_global_atoms(mtop);
        t_block mols  = gmx_mtop_molecules_t_block(mtop);

        tpx_make_chain_identifiers(&atoms, &mols);
        checkChainIds(atoms, expectedChains);

        done_atom(&atoms);
        done_block(&mols);
    }

    //! Test tpx_make_chain_identifiers for input legacy top
    static void runTopologyTest(t_topology& top, ArrayRef<const ParticleChain> expectedChains)
    {
        tpx_make_chain_identifiers(&top.atoms, &top.mols);
        checkChainIds(top.atoms, expectedChains);
    }

    //! Check chain IDs in \p atoms against reference data
    static void checkChainIds(const t_atoms& atoms, ArrayRef<const ParticleChain> expectedChains)
    {
        std::vector<ParticleChain> chains{ { 1, atoms.resinfo[atoms.atom[0].resind].chainid } };
        for (int i = 1; i < atoms.nr; ++i)
        {
            if (const char nextChainId = atoms.resinfo[atoms.atom[i].resind].chainid;
                nextChainId != chains.back().id)
            {
                chains.push_back({ 0, nextChainId });
            }
            ++chains.back().length;
        }

        EXPECT_EQ(chains.size(), expectedChains.size());
        for (int i = 0; i < gmx::ssize(expectedChains); ++i)
        {
            EXPECT_EQ(chains[i].id, expectedChains[i].id);
            EXPECT_EQ(chains[i].length, expectedChains[i].length);
        }
    }
};

TEST_F(ChainIdFromTopologyTest, WithFullTopology)
{
    // OctaneSandwich fulfils current chain ID construction criteria:
    // >=2 mols with >=15 atoms each
    TprAndFileManager tprManager("OctaneSandwich");
    tprManager.setInputConf("OctaneSandwich.gro");
    tprManager.setPositionRestraintConf("OctaneSandwich.gro");
    runTest(tprManager.tprName(),
            std::vector<ParticleChain>{ { 166, 'A' }, { 798, ' ' }, { 166, 'B' }, { 12691, ' ' } });
}

TEST_F(ChainIdFromTopologyTest, WithoutFullTopology)
{
    // OctaneSandwich fulfils current chain ID construction criteria:
    // >=2 mols with >=15 atoms each
    TestFileManager fileManager;
    runTest(fileManager.getInputFilePath("OctaneSandwich.gro"),
            std::vector<ParticleChain>{ { 13821, ' ' } });
}

TEST_F(ChainIdFromTopologyTest, NoChainsForSingleMoleculeTopology)
{
    // lysozyme does not fulfil current chain ID construction criteria:
    // >=2 mols with >=15 atoms each (only has 1 molecule with >= 15 atoms)
    TprAndFileManager tprManager("lysozyme");
    runTest(tprManager.tprName(), std::vector<ParticleChain>{ { 156, ' ' } });
}

} // namespace
} // namespace test
} // namespace gmx
