/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * Implements test of index generation routines
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_topology
 */
#include "gmxpre.h"

#include "gromacs/topology/index.h"

#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/cmdlinetest.h"
#include "testutils/filematchers.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

void checkBlocks(TestReferenceChecker* checker, const t_blocka* blocks)
{
    TestReferenceChecker compound(checker->checkCompound("Blocks", nullptr));
    compound.checkInteger(blocks->nr, "Number");
    compound.checkSequence(blocks->index, blocks->index + blocks->nr, "Index");
    compound.checkInteger(blocks->nra, "NumberAtoms");
    compound.checkSequence(blocks->a, blocks->a + blocks->nra, "AtomIndex");
}

void compareBlocks(const t_blocka* one, const t_blocka* two)
{
    ASSERT_EQ(one->nr, two->nr);
    for (int i = 0; i < one->nr; ++i)
    {
        EXPECT_EQ(one->index[i], two->index[i]);
    }
    ASSERT_EQ(one->nra, two->nra);
    for (int i = 0; i < one->nra; ++i)
    {
        EXPECT_EQ(one->a[i], two->a[i]);
    }
}

void checkFileMatch(TestReferenceChecker* checker, const std::string& fileName, const std::string& fullPath)
{
    TestReferenceChecker fileChecker(checker->checkCompound("File", fileName.c_str()));
    auto                 matcher = TextFileMatch(ExactTextMatch()).createFileMatcher();
    matcher->checkFile(fullPath, &fileChecker);
}

class IndexTest : public ::testing::Test
{
public:
    IndexTest();
    ~IndexTest() override;

    //! Handle to atoms from topology.
    const t_atoms* atoms() { return topInfo_.atoms(); }
    //! Handle to checker.
    TestReferenceChecker* checker() { return &checker_; }
    //! Handle to file manager.
    TestFileManager* manager() { return &manager_; }
    //! Index group structure.
    t_blocka blocks_;
    //! Index group names.
    char** groupNames_ = nullptr;

private:
    //! Input structure data.
    TopologyInformation topInfo_;
    //! File manager for test.
    TestFileManager manager_;
    //! Handler for reference data.
    TestReferenceData data_;
    //! Handler for checking test data.
    TestReferenceChecker checker_;
};

IndexTest::IndexTest() : checker_(data_.rootChecker())
{
    // When we have many test cases using this class, refactor to fill
    // a static topInfo only once, in SetUpTestSuite()
    topInfo_.fillFromInputFile(manager()->getInputFilePath("lysozyme.gro"));
    init_blocka(&blocks_);
    snew(groupNames_, 1);
}

IndexTest::~IndexTest()
{
    for (int i = 0; i < blocks_.nr; ++i)
    {
        sfree(groupNames_[i]);
    }
    sfree(groupNames_);
    done_blocka(&blocks_);
}

TEST_F(IndexTest, AnalyseWorksDefaultGroups)
{
    analyse(atoms(), &blocks_, &groupNames_, false, false);
    checkBlocks(checker(), &blocks_);
}

TEST_F(IndexTest, WriteIndexWorks)
{
    analyse(atoms(), &blocks_, &groupNames_, false, false);
    std::string fileName = "out.ndx";
    std::string fullPath = manager()->getTemporaryFilePath(fileName);
    write_index(fullPath.c_str(), &blocks_, groupNames_, false, atoms()->nr);
    checkFileMatch(checker(), fileName, fullPath);
}

TEST_F(IndexTest, WriteAndReadIndexWorks)
{
    analyse(atoms(), &blocks_, &groupNames_, false, false);
    std::string fileName = "out.ndx";
    std::string fullPath = manager()->getTemporaryFilePath(fileName);
    write_index(fullPath.c_str(), &blocks_, groupNames_, false, atoms()->nr);
    char**    newNames = nullptr;
    t_blocka* newIndex = init_index(fullPath.c_str(), &newNames);
    compareBlocks(&blocks_, newIndex);
    for (int i = 0; i < newIndex->nr; ++i)
    {
        sfree(newNames[i]);
    }
    sfree(newNames);
    done_blocka(newIndex);
    sfree(newIndex);
}

} // namespace

} // namespace test

} // namespace gmx
