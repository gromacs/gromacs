/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * Tests for TopologyInformation
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/topologyinformation.h"

#include <filesystem>
#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"


namespace gmx
{
namespace test
{
namespace
{

TEST(TopologyInformation, CantWorkWithoutReadingAFile)
{
    TopologyInformation topInfo;
    EXPECT_FALSE(topInfo.hasTopology());
    EXPECT_FALSE(topInfo.hasFullTopology());
    EXPECT_EQ(nullptr, topInfo.mtop());
    EXPECT_EQ(nullptr, topInfo.expandedTopology());
    auto atoms1 = topInfo.copyAtoms();
    EXPECT_TRUE(atoms1);
    auto atoms2 = topInfo.copyAtoms();
    ASSERT_TRUE(atoms2);
    EXPECT_NE(atoms1.get(), atoms2.get());
    EXPECT_EQ(0, atoms1->nr);
    EXPECT_EQ(PbcType::Unset, topInfo.pbcType());
    EXPECT_THROW(topInfo.x().size(), gmx::APIError);
    EXPECT_THROW(topInfo.v().size(), gmx::APIError);
    matrix box{ { -2 } };
    topInfo.getBox(box);
    EXPECT_EQ(0, box[XX][XX]);
    EXPECT_EQ(0, box[XX][YY]);
    EXPECT_EQ(0, box[XX][ZZ]);
    EXPECT_EQ(0, box[YY][XX]);
    EXPECT_EQ(0, box[YY][YY]);
    EXPECT_EQ(0, box[YY][ZZ]);
    EXPECT_EQ(0, box[ZZ][XX]);
    EXPECT_EQ(0, box[ZZ][YY]);
    EXPECT_EQ(0, box[ZZ][ZZ]);
    EXPECT_FALSE(topInfo.name());
}

//! Common test code to reduce duplication
void runCommonTests(const TopologyInformation& topInfo, const int numAtoms)
{
    EXPECT_TRUE(topInfo.hasTopology());
    ASSERT_TRUE(topInfo.mtop());
    EXPECT_EQ(numAtoms, topInfo.mtop()->natoms);
    // TODO Dump mtop to refdata when that is possible
    ASSERT_TRUE(topInfo.expandedTopology());
    auto atoms1 = topInfo.copyAtoms();
    EXPECT_TRUE(atoms1);
    auto atoms2 = topInfo.copyAtoms();
    EXPECT_TRUE(atoms2);
    // Must be different pointer to a deep copy.
    EXPECT_NE(atoms1.get(), atoms2.get());
    const auto* atoms = topInfo.atoms();
    // Must be a pointer to a different instance.
    EXPECT_NE(atoms1.get(), atoms);
    EXPECT_NE(atoms2.get(), atoms);
    EXPECT_EQ(numAtoms, topInfo.x().size());
    EXPECT_EQ(numAtoms, topInfo.v().size());
    matrix box{ { -2 } };
    topInfo.getBox(box);
    EXPECT_FLOAT_EQ(5.9062, box[XX][XX]);
    EXPECT_FLOAT_EQ(0, box[XX][YY]);
    EXPECT_FLOAT_EQ(0, box[XX][ZZ]);
    EXPECT_FLOAT_EQ(0, box[YY][XX]);
    EXPECT_FLOAT_EQ(6.8451, box[YY][YY]);
    EXPECT_FLOAT_EQ(0, box[YY][ZZ]);
    EXPECT_FLOAT_EQ(0, box[ZZ][XX]);
    EXPECT_FLOAT_EQ(0, box[ZZ][YY]);
    EXPECT_FLOAT_EQ(3.0517, box[ZZ][ZZ]);
    EXPECT_STREQ("First 10 residues from 1AKI", topInfo.name());
}

TEST(TopologyInformation, WorksWithGroFile)
{
    const int           numAtoms = 156;
    TopologyInformation topInfo;
    topInfo.fillFromInputFile(TestFileManager::getInputFilePath("lysozyme.gro").string());
    EXPECT_FALSE(topInfo.hasFullTopology());
    runCommonTests(topInfo, numAtoms);
    EXPECT_EQ(PbcType::Unset, topInfo.pbcType());

    // Check the per-atom data
    auto atoms = topInfo.copyAtoms();
    ASSERT_EQ(numAtoms, atoms->nr);
    EXPECT_TRUE(atoms->haveMass);
    // TODO atommass.dat assumes united atom CA, which is probably not expected behaviour
    EXPECT_FLOAT_EQ(13.019, atoms->atom[26].m);
    EXPECT_FALSE(atoms->haveCharge);
    EXPECT_FALSE(atoms->haveType);
    EXPECT_EQ(0, atoms->atom[26].type);
    EXPECT_EQ(0, atoms->atom[26].atomnumber);
    EXPECT_EQ(1, atoms->atom[26].resind);
    // gro files don't have the information that pdb files might
    EXPECT_FALSE(atoms->havePdbInfo);
    EXPECT_FALSE(atoms->pdbinfo);
    EXPECT_EQ(10, atoms->nres);
    ASSERT_TRUE(atoms->resinfo);
    ASSERT_TRUE(atoms->resinfo[4].name);
    EXPECT_STREQ("ARG", *atoms->resinfo[4].name);
    EXPECT_EQ(5, atoms->resinfo[4].nr);
    EXPECT_EQ(0, atoms->resinfo[4].chainnum);
    EXPECT_EQ(' ', atoms->resinfo[4].chainid);
}

TEST(TopologyInformation, WorksWithPdbFile)
{
    const int           numAtoms = 156;
    TopologyInformation topInfo;
    topInfo.fillFromInputFile(TestFileManager::getInputFilePath("lysozyme.pdb").string());
    EXPECT_FALSE(topInfo.hasFullTopology());
    runCommonTests(topInfo, numAtoms);
    // TODO why does this differ from .gro?
    EXPECT_EQ(PbcType::Xyz, topInfo.pbcType());

    // Check the per-atom data
    auto atoms = topInfo.copyAtoms();
    ASSERT_EQ(numAtoms, atoms->nr);
    EXPECT_TRUE(atoms->haveMass);
    // TODO atommass.dat assumes united atom CA, which is probably not expected behaviour
    EXPECT_FLOAT_EQ(13.019, atoms->atom[26].m);
    EXPECT_FALSE(atoms->haveCharge);
    EXPECT_FALSE(atoms->haveType);
    EXPECT_EQ(0, atoms->atom[26].type);
    EXPECT_EQ(0, atoms->atom[26].atomnumber);
    EXPECT_EQ(1, atoms->atom[26].resind);
    // pdb files can carry more information than gro
    EXPECT_TRUE(atoms->havePdbInfo);
    ASSERT_TRUE(atoms->pdbinfo);
    EXPECT_EQ(10, atoms->nres);
    ASSERT_TRUE(atoms->resinfo);
    ASSERT_TRUE(atoms->resinfo[4].name);
    EXPECT_STREQ("ARG", *atoms->resinfo[4].name);
    EXPECT_EQ(5, atoms->resinfo[4].nr);
    EXPECT_EQ(0, atoms->resinfo[4].chainnum);
    EXPECT_EQ('B', atoms->resinfo[4].chainid);
}

TEST(TopologyInformation, WorksWithTprFromPdbFile)
{
    TestFileManager fileManager;

    // Make the tpr file to use
    std::string       name             = "lysozyme";
    const std::string mdpInputFileName = fileManager.getTemporaryFilePath(name + ".mdp").string();
    // Ensure the seeds have a value so that the resulting .tpr dump
    // is reproducible.
    TextWriter::writeFileFromString(mdpInputFileName, "");
    std::string tprName = fileManager.getTemporaryFilePath(name + ".tpr").string();
    {
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-f", mdpInputFileName);
        caller.addOption("-p", TestFileManager::getInputFilePath(name + ".top").string());
        caller.addOption("-c", TestFileManager::getInputFilePath(name + ".pdb").string());
        caller.addOption("-o", tprName);
        ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
    }

    const int           numAtoms = 156;
    TopologyInformation topInfo;
    topInfo.fillFromInputFile(tprName);
    EXPECT_TRUE(topInfo.hasFullTopology());
    runCommonTests(topInfo, numAtoms);
    // TODO why does this differ from .gro?
    EXPECT_EQ(PbcType::Xyz, topInfo.pbcType());

    // Check the per-atom data
    auto atoms = topInfo.copyAtoms();
    ASSERT_EQ(numAtoms, atoms->nr);
    EXPECT_TRUE(atoms->haveMass);
    EXPECT_FLOAT_EQ(12.011, atoms->atom[26].m);
    EXPECT_TRUE(atoms->haveCharge);
    EXPECT_TRUE(atoms->haveType);
    EXPECT_EQ(2, atoms->atom[26].type);
    EXPECT_EQ(6, atoms->atom[26].atomnumber);
    EXPECT_EQ(1, atoms->atom[26].resind);
    // tpr files also don't carry pdb information
    EXPECT_FALSE(atoms->havePdbInfo);
    EXPECT_FALSE(atoms->pdbinfo);
    EXPECT_EQ(10, atoms->nres);
    ASSERT_TRUE(atoms->resinfo);
    ASSERT_TRUE(atoms->resinfo[4].name);
    EXPECT_STREQ("ARG", *atoms->resinfo[4].name);
    EXPECT_EQ(5, atoms->resinfo[4].nr);
    EXPECT_EQ(0, atoms->resinfo[4].chainnum);
    // In particular, chain ID does not get recorded in the .tpr file
    EXPECT_EQ(0, atoms->resinfo[4].chainid);
}

} // namespace
} // namespace test
} // namespace gmx
