/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
 * Tests the index group handling in the selection engine.
 *
 * \todo
 * Tests for other functions, at least the set operations.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/selection/indexutil.h"

#include <gtest/gtest.h>

#include "gromacs/topology/block.h"

#include "testutils/refdata.h"

#include "toputils.h"

namespace
{

/********************************************************************
 * IndexBlockTest
 */

class IndexBlockTest : public ::testing::Test
{
    public:
        IndexBlockTest();
        ~IndexBlockTest();

        void setGroup(int count, const int atoms[]);
        template <int count>
        void setGroup(const int (&atoms)[count])
        {
            setGroup(count, atoms);
        }

        void checkBlocka();

        gmx::test::TestReferenceData    data_;
        gmx::test::TopologyManager      topManager_;
        gmx_ana_index_t                 g_;
        t_blocka                        blocka_;
};

IndexBlockTest::IndexBlockTest()
{
    blocka_.nr           = 0;
    blocka_.index        = NULL;
    blocka_.nalloc_index = 0;
    blocka_.nra          = 0;
    blocka_.a            = NULL;
    blocka_.nalloc_a     = 0;
    gmx_ana_index_clear(&g_);
}

IndexBlockTest::~IndexBlockTest()
{
    done_blocka(&blocka_);
}

void IndexBlockTest::setGroup(int count, const int atoms[])
{
    g_.isize = count;
    g_.index = const_cast<int *>(atoms);
}

void IndexBlockTest::checkBlocka()
{
    gmx::test::TestReferenceChecker compound(
            data_.rootChecker().checkCompound("BlockAtoms", "Block"));
    compound.checkSequenceArray(g_.isize, g_.index, "Input");
    compound.checkInteger(blocka_.nr, "Count");
    for (int i = 0; i < blocka_.nr; ++i)
    {
        gmx::test::TestReferenceChecker blockCompound(
                compound.checkCompound("Block", NULL));
        blockCompound.checkSequence(&blocka_.a[blocka_.index[i]],
                                    &blocka_.a[blocka_.index[i+1]],
                                    "Atoms");
    }
}

/********************************************************************
 * gmx_ana_index_make_block() tests
 */

TEST_F(IndexBlockTest, CreatesUnknownBlock)
{
    gmx_ana_index_make_block(&blocka_, NULL, NULL, INDEX_UNKNOWN, false);
    checkBlocka();
    done_blocka(&blocka_);
    gmx_ana_index_make_block(&blocka_, NULL, NULL, INDEX_UNKNOWN, false);
    checkBlocka();
}

TEST_F(IndexBlockTest, CreatesAtomBlock)
{
    const int group[] = { 0, 1, 3, 4, 6 };
    setGroup(group);
    gmx_ana_index_make_block(&blocka_, NULL, &g_, INDEX_ATOM, false);
    checkBlocka();
    done_blocka(&blocka_);
    gmx_ana_index_make_block(&blocka_, NULL, &g_, INDEX_ATOM, true);
    checkBlocka();
}

TEST_F(IndexBlockTest, CreatesResidueBlock)
{
    const int group[] = { 0, 1, 3, 6, 7 };
    topManager_.initAtoms(9);
    topManager_.initUniformResidues(3);
    setGroup(group);
    gmx_ana_index_make_block(&blocka_, topManager_.topology(), &g_,
                             INDEX_RES, false);
    checkBlocka();
}

TEST_F(IndexBlockTest, CreatesMoleculeBlock)
{
    const int group[] = { 3, 4, 7, 8, 13 };
    topManager_.initAtoms(18);
    topManager_.initUniformMolecules(3);
    setGroup(group);
    gmx_ana_index_make_block(&blocka_, topManager_.topology(), &g_,
                             INDEX_MOL, false);
    checkBlocka();
}

TEST_F(IndexBlockTest, CreatesResidueBlockWithCompletion)
{
    const int group[] = { 3, 4, 7, 8, 13 };
    topManager_.initAtoms(18);
    topManager_.initUniformResidues(3);
    setGroup(group);
    gmx_ana_index_make_block(&blocka_, topManager_.topology(), &g_,
                             INDEX_RES, true);
    checkBlocka();
}

TEST_F(IndexBlockTest, CreatesMoleculeBlockWithCompletion)
{
    const int group[] = { 3, 4, 7, 8, 13 };
    topManager_.initAtoms(18);
    topManager_.initUniformMolecules(3);
    setGroup(group);
    gmx_ana_index_make_block(&blocka_, topManager_.topology(), &g_,
                             INDEX_MOL, true);
    checkBlocka();
}

TEST_F(IndexBlockTest, CreatesSingleBlock)
{
    const int group[] = { 0, 1, 3, 4, 6 };
    setGroup(group);
    gmx_ana_index_make_block(&blocka_, NULL, &g_, INDEX_ALL, false);
    checkBlocka();
    done_blocka(&blocka_);
    gmx_ana_index_make_block(&blocka_, NULL, &g_, INDEX_ALL, true);
    checkBlocka();
}

/********************************************************************
 * gmx_ana_index_has_full_ablocks() tests
 */

TEST_F(IndexBlockTest, ChecksGroupForFullBlocksPositive)
{
    const int maxGroup[] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17
    };
    const int testGroup[] = { 3, 4, 5, 6, 7, 8, 12, 13, 14 };
    topManager_.initAtoms(18);
    topManager_.initUniformResidues(3);
    setGroup(maxGroup);
    gmx_ana_index_make_block(&blocka_, topManager_.topology(), &g_,
                             INDEX_RES, false);
    setGroup(testGroup);
    EXPECT_TRUE(gmx_ana_index_has_full_ablocks(&g_, &blocka_));
}

TEST_F(IndexBlockTest, ChecksOutOfOrderGroupForFullBlocksPositive)
{
    const int maxGroup[] = {
        15, 16, 17, 2, 1, 0, 12, 13, 14, 5, 4, 3, 9, 10, 11, 8, 7, 6
    };
    const int testGroup[] = { 2, 1, 0, 5, 4, 3, 8, 7, 6, };
    topManager_.initAtoms(18);
    topManager_.initUniformResidues(3);
    setGroup(maxGroup);
    gmx_ana_index_make_block(&blocka_, topManager_.topology(), &g_,
                             INDEX_RES, false);
    setGroup(testGroup);
    EXPECT_TRUE(gmx_ana_index_has_full_ablocks(&g_, &blocka_));
}

TEST_F(IndexBlockTest, ChecksGroupForFullBlocksNegative)
{
    const int maxGroup[] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17
    };
    const int testGroup1[] = { 3, 4, 5, 6, 7, 8, 12, 13 };
    const int testGroup2[] = { 3, 4, 5, 6, 7, 12, 13, 14 };
    const int testGroup3[] = { 4, 5, 6, 7, 8, 12, 13, 14 };

    topManager_.initAtoms(18);
    topManager_.initUniformResidues(3);
    setGroup(maxGroup);
    gmx_ana_index_make_block(&blocka_, topManager_.topology(), &g_,
                             INDEX_RES, false);

    setGroup(testGroup1);
    EXPECT_FALSE(gmx_ana_index_has_full_ablocks(&g_, &blocka_));

    setGroup(testGroup2);
    EXPECT_FALSE(gmx_ana_index_has_full_ablocks(&g_, &blocka_));

    setGroup(testGroup3);
    EXPECT_FALSE(gmx_ana_index_has_full_ablocks(&g_, &blocka_));
}

/********************************************************************
 * gmx_ana_index_has_complete_elems() tests
 */

TEST_F(IndexBlockTest, ChecksGroupForCompleteElementsTrivial)
{
    const int group[] = { 0, 1, 2 };
    setGroup(group);
    EXPECT_TRUE(gmx_ana_index_has_complete_elems(&g_, INDEX_ATOM, NULL));
    EXPECT_FALSE(gmx_ana_index_has_complete_elems(&g_, INDEX_ALL, NULL));
    EXPECT_FALSE(gmx_ana_index_has_complete_elems(&g_, INDEX_UNKNOWN, NULL));
}

TEST_F(IndexBlockTest, ChecksGroupForCompleteResiduesPositive)
{
    const int group1[] = { 0, 1, 2, 6, 7, 8, 12, 13, 14 };
    const int group2[] = { 3, 4, 5, 6, 7, 8 };

    topManager_.initAtoms(15);
    topManager_.initUniformResidues(3);
    t_topology *top = topManager_.topology();

    setGroup(group1);
    EXPECT_TRUE(gmx_ana_index_has_complete_elems(&g_, INDEX_RES, top));

    setGroup(group2);
    EXPECT_TRUE(gmx_ana_index_has_complete_elems(&g_, INDEX_RES, top));
}

TEST_F(IndexBlockTest, ChecksGroupForCompleteResiduesNegative)
{
    const int group1[] = { 3, 4, 5, 6, 7, 8, 12, 13 };
    const int group2[] = { 3, 4, 5, 6, 7, 12, 13, 14 };
    const int group3[] = { 4, 5, 6, 7, 8, 12, 13, 14 };

    topManager_.initAtoms(18);
    topManager_.initUniformResidues(3);
    t_topology *top = topManager_.topology();

    setGroup(group1);
    EXPECT_FALSE(gmx_ana_index_has_complete_elems(&g_, INDEX_RES, top));

    setGroup(group2);
    EXPECT_FALSE(gmx_ana_index_has_complete_elems(&g_, INDEX_RES, top));

    setGroup(group3);
    EXPECT_FALSE(gmx_ana_index_has_complete_elems(&g_, INDEX_RES, top));
}

TEST_F(IndexBlockTest, ChecksGroupForCompleteMoleculesPositive)
{
    const int group[] = { 0, 1, 2, 6, 7, 8, 12, 13, 14 };

    topManager_.initAtoms(15);
    topManager_.initUniformMolecules(3);
    t_topology *top = topManager_.topology();

    setGroup(group);
    EXPECT_TRUE(gmx_ana_index_has_complete_elems(&g_, INDEX_MOL, top));
}

TEST_F(IndexBlockTest, ChecksGroupForCompleteMoleculesNegative)
{
    const int group1[] = { 3, 4, 5, 6, 7, 8, 12, 13 };
    const int group2[] = { 3, 4, 5, 6, 7, 12, 13, 14 };
    const int group3[] = { 4, 5, 6, 7, 8, 12, 13, 14 };

    topManager_.initAtoms(18);
    topManager_.initUniformMolecules(3);
    t_topology *top = topManager_.topology();

    setGroup(group1);
    EXPECT_FALSE(gmx_ana_index_has_complete_elems(&g_, INDEX_MOL, top));

    setGroup(group2);
    EXPECT_FALSE(gmx_ana_index_has_complete_elems(&g_, INDEX_MOL, top));

    setGroup(group3);
    EXPECT_FALSE(gmx_ana_index_has_complete_elems(&g_, INDEX_MOL, top));
}

/********************************************************************
 * IndexMapTest
 */

class IndexMapTest : public ::testing::Test
{
    public:
        IndexMapTest();
        ~IndexMapTest();

        void testInit(int atomCount, const int atoms[], e_index_t type);
        void testUpdate(int atomCount, const int atoms[], bool bMaskOnly,
                        const char *name);
        void testOrgIdGroup(e_index_t type, const char *name);
        template <int count>
        void testInit(const int (&atoms)[count], e_index_t type)
        {
            testInit(count, atoms, type);
        }
        template <int count>
        void testUpdate(const int (&atoms)[count], bool bMaskOnly,
                        const char *name)
        {
            testUpdate(count, atoms, bMaskOnly, name);
        }

        void checkMapping(int atomCount, const int atoms[], const char *name);

        gmx::test::TestReferenceData    data_;
        gmx::test::TestReferenceChecker checker_;
        gmx::test::TopologyManager      topManager_;
        gmx_ana_indexmap_t              map_;

    private:
        gmx_ana_index_t                 initGroup_;
};

IndexMapTest::IndexMapTest()
    : checker_(data_.rootChecker())
{
    gmx_ana_indexmap_clear(&map_);
    gmx_ana_index_clear(&initGroup_);
}

IndexMapTest::~IndexMapTest()
{
    gmx_ana_indexmap_deinit(&map_);
}

void IndexMapTest::testInit(int atomCount, const int atoms[], e_index_t type)
{
    initGroup_.isize = atomCount;
    initGroup_.index = const_cast<int *>(atoms);
    gmx_ana_indexmap_init(&map_, &initGroup_, topManager_.topology(), type);
    EXPECT_EQ(type, map_.type);
    checkMapping(atomCount, atoms, "Initialized");
}

void IndexMapTest::testUpdate(int atomCount, const int atoms[], bool bMaskOnly,
                              const char *name)
{
    gmx_ana_index_t g;
    g.isize = atomCount;
    g.index = const_cast<int *>(atoms);
    gmx_ana_indexmap_update(&map_, &g, bMaskOnly);
    if (name == NULL)
    {
        name = "Updated";
    }
    if (bMaskOnly)
    {
        checkMapping(initGroup_.isize, initGroup_.index, name);
    }
    else
    {
        checkMapping(atomCount, atoms, name);
    }
}

void IndexMapTest::testOrgIdGroup(e_index_t type, const char *name)
{
    gmx::test::TestReferenceChecker compound(
            checker_.checkCompound("OrgIdGroups", name));
    const int count
        = gmx_ana_indexmap_init_orgid_group(&map_, topManager_.topology(), type);
    compound.checkInteger(count, "GroupCount");
    compound.checkSequenceArray(map_.mapb.nr, map_.orgid, "OrgId");
    for (int i = 0; i < map_.mapb.nr; ++i)
    {
        EXPECT_EQ(map_.orgid[i], map_.mapid[i]);
    }
}

void IndexMapTest::checkMapping(int atomCount, const int atoms[],
                                const char *name)
{
    gmx::test::TestReferenceChecker compound(
            checker_.checkCompound("IndexMapping", name));
    compound.checkSequenceArray(atomCount, atoms, "Input");
    compound.checkInteger(map_.mapb.nr, "Count");
    for (int i = 0; i < map_.mapb.nr; ++i)
    {
        gmx::test::TestReferenceChecker blockCompound(
                compound.checkCompound("Block", NULL));
        blockCompound.checkSequence(&atoms[map_.mapb.index[i]],
                                    &atoms[map_.mapb.index[i+1]],
                                    "Atoms");
        blockCompound.checkInteger(map_.refid[i], "RefId");
        blockCompound.checkInteger(map_.mapid[i], "MapId");
        int originalIdIndex = (map_.refid[i] != -1 ? map_.refid[i] : i);
        EXPECT_EQ(map_.orgid[originalIdIndex], map_.mapid[i]);
    }
}

/********************************************************************
 * gmx_ana_indexmap_t tests
 */

TEST_F(IndexMapTest, InitializesAtomBlock)
{
    const int maxGroup[] = { 1, 2, 4, 5 };
    testInit(maxGroup, INDEX_ATOM);
}

TEST_F(IndexMapTest, InitializesOrgIdGroupAtom)
{
    const int maxGroup[] = { 2, 5, 7 };
    testInit(maxGroup, INDEX_ATOM);
    testOrgIdGroup(INDEX_ATOM, "Atoms");
}

TEST_F(IndexMapTest, InitializesOrgIdGroupSingle)
{
    const int maxGroup[] = { 3, 4, 7, 8, 13 };
    topManager_.initAtoms(18);
    topManager_.initUniformResidues(3);
    testInit(maxGroup, INDEX_RES);
    testOrgIdGroup(INDEX_ATOM, "Single");
}

TEST_F(IndexMapTest, InitializesOrgIdGroupResidue)
{
    const int maxGroup[] = { 3, 4, 7, 8, 13 };
    topManager_.initAtoms(18);
    topManager_.initUniformResidues(3);
    testInit(maxGroup, INDEX_ATOM);
    testOrgIdGroup(INDEX_RES, "Residues");
}

TEST_F(IndexMapTest, InitializesOrgIdGroupMolecule)
{
    const int maxGroup[] = { 1, 2, 3, 4, 7, 8, 13 };
    topManager_.initAtoms(18);
    topManager_.initUniformResidues(3);
    topManager_.initUniformMolecules(6);
    testInit(maxGroup, INDEX_RES);
    testOrgIdGroup(INDEX_MOL, "Molecules");
}

TEST_F(IndexMapTest, InitializesOrgIdGroupAll)
{
    const int maxGroup[] = { 3, 4, 7, 8, 13 };
    testInit(maxGroup, INDEX_ATOM);
    testOrgIdGroup(INDEX_ALL, "All");
}

TEST_F(IndexMapTest, InitializesMoleculeBlock)
{
    const int maxGroup[] = { 3, 4, 7, 8, 13 };
    topManager_.initAtoms(18);
    topManager_.initUniformMolecules(3);
    testInit(maxGroup, INDEX_MOL);
}

TEST_F(IndexMapTest, MapsSingleBlock)
{
    const int maxGroup[]  = { 0, 1, 2, 3 };
    const int evalGroup[] = { 0, 2 };
    testInit(maxGroup, INDEX_ALL);
    testUpdate(evalGroup, false, NULL);
}

TEST_F(IndexMapTest, MapsResidueBlocks)
{
    const int maxGroup[] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17
    };
    const int evalGroup[] = { 3, 4, 7, 8, 13 };
    topManager_.initAtoms(18);
    topManager_.initUniformResidues(3);
    testInit(maxGroup, INDEX_RES);
    testUpdate(evalGroup, false, NULL);
}

TEST_F(IndexMapTest, MapsResidueBlocksWithMask)
{
    const int maxGroup[] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17
    };
    const int evalGroup[] = { 3, 4, 7, 8, 13 };
    topManager_.initAtoms(18);
    topManager_.initUniformResidues(3);
    testInit(maxGroup, INDEX_RES);
    testUpdate(evalGroup, true, NULL);
}

TEST_F(IndexMapTest, HandlesMultipleRequests)
{
    const int maxGroup[] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17
    };
    const int evalGroup[] = { 3, 4, 7, 8, 13 };
    topManager_.initAtoms(18);
    topManager_.initUniformResidues(3);
    testInit(maxGroup, INDEX_RES);
    testUpdate(evalGroup, false, "EvaluatedNoMask");
    testUpdate(evalGroup, true, "EvaluatedMask");
    testUpdate(maxGroup, true, "Initialized");
    testUpdate(evalGroup, true, "EvaluatedMask");
    testUpdate(evalGroup, false, "EvaluatedNoMask");
    testUpdate(maxGroup, false, "Initialized");
}

} // namespace
