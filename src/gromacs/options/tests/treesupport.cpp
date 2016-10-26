/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * Tests option support for operations on KeyValueTree.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "gromacs/options/treesupport.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/options/repeatingsection.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace
{

/********************************************************************
 * Tests for assignOptionsFromKeyValueTree()
 */

TEST(TreeValueSupportAssignTest, AssignsFromTree)
{
    int                      a0 = 0, a1 = 0;
    std::string              b1;

    gmx::Options             options;
    options.addOption(gmx::IntegerOption("a").store(&a0));
    auto                     sec = options.addSection(gmx::OptionSection("s"));
    sec.addOption(gmx::IntegerOption("a").store(&a1));
    sec.addOption(gmx::StringOption("b").store(&b1));

    gmx::KeyValueTreeBuilder builder;
    builder.rootObject().addValue<int>("a", 2);
    auto                     obj = builder.rootObject().addObject("s");
    obj.addValue<int>("a", 1);
    obj.addValue<std::string>("b", "foo");
    gmx::KeyValueTreeObject  tree = builder.build();

    ASSERT_NO_THROW_GMX(gmx::assignOptionsFromKeyValueTree(&options, tree, nullptr));
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_EQ(2, a0);
    EXPECT_EQ(1, a1);
    EXPECT_EQ("foo", b1);
}

struct SectionData
{
    int a;
};

TEST(TreeValueSupportAssignTest, AssignsFromTreeWithArrays)
{
    std::vector<int>         a0;
    std::vector<SectionData> s;

    gmx::Options             options;
    options.addOption(gmx::IntegerOption("a").storeVector(&a0).multiValue());
    auto                     sec = options.addSection(gmx::RepeatingOptionSection<SectionData>("s").storeVector(&s));
    sec.addOption(gmx::IntegerOption("a").store(&sec.bind().a));

    gmx::KeyValueTreeBuilder builder;
    auto                     array = builder.rootObject().addUniformArray<int>("a");
    array.addValue(1);
    array.addValue(2);
    auto                     objArray = builder.rootObject().addObjectArray("s");
    auto                     obj1     = objArray.addObject();
    obj1.addValue<int>("a", 3);
    auto                     obj2 = objArray.addObject();
    obj2.addValue<int>("a", 4);
    gmx::KeyValueTreeObject  tree = builder.build();

    ASSERT_NO_THROW_GMX(gmx::assignOptionsFromKeyValueTree(&options, tree, nullptr));
    EXPECT_NO_THROW_GMX(options.finish());

    ASSERT_EQ(2U, a0.size());
    EXPECT_EQ(1, a0[0]);
    EXPECT_EQ(2, a0[1]);
    ASSERT_EQ(2U, s.size());
    EXPECT_EQ(3, s[0].a);
    EXPECT_EQ(4, s[1].a);
}

TEST(TreeValueSupportAssignErrorTest, HandlesInvalidValue)
{
    int                      a1 = 0;

    gmx::Options             options;
    auto                     sec = options.addSection(gmx::OptionSection("s"));
    sec.addOption(gmx::IntegerOption("a").store(&a1));

    gmx::KeyValueTreeBuilder builder;
    auto                     obj = builder.rootObject().addObject("s");
    obj.addValue<std::string>("a", "foo");
    gmx::KeyValueTreeObject  tree = builder.build();

    EXPECT_THROW_GMX(gmx::assignOptionsFromKeyValueTree(&options, tree, nullptr),
                     gmx::InvalidInputError);
}

/********************************************************************
 * Tests for adjustKeyValueTreeFromOptions()
 */

class TreeValueSupportAdjustTest : public ::testing::Test
{
    public:
        void runTest()
        {
            gmx::test::TestReferenceData    refdata;
            gmx::test::TestReferenceChecker checker(refdata.rootChecker());
            gmx::KeyValueTreeObject         tree(builder_.build());
            checker.checkKeyValueTreeObject(tree, "Input");
            ASSERT_NO_THROW_GMX(tree = gmx::adjustKeyValueTreeFromOptions(tree, options_));
            checker.checkKeyValueTreeObject(tree, "Output");
        }

        gmx::Options              options_;
        gmx::KeyValueTreeBuilder  builder_;
};

TEST_F(TreeValueSupportAdjustTest, FillsDefaultValues)
{
    options_.addOption(gmx::IntegerOption("a").defaultValue(2));
    runTest();
}

TEST_F(TreeValueSupportAdjustTest, FillsDefaultVectorValues)
{
    int v[3] = {1, 2, 3};
    options_.addOption(gmx::IntegerOption("a").store(v).vector());
    runTest();
}

TEST_F(TreeValueSupportAdjustTest, FillsDefaultObjectValues)
{
    auto sec1 = options_.addSection(gmx::OptionSection("s"));
    sec1.addOption(gmx::IntegerOption("a").defaultValue(1));
    auto sec2 = options_.addSection(gmx::OptionSection("r"));
    sec2.addOption(gmx::IntegerOption("a").defaultValue(2));
    options_.addOption(gmx::IntegerOption("a").defaultValue(3));
    runTest();
}

TEST_F(TreeValueSupportAdjustTest, NormalizesValues)
{
    options_.addOption(gmx::IntegerOption("a"));
    builder_.rootObject().addValue<std::string>("a", "2");
    runTest();
}

TEST_F(TreeValueSupportAdjustTest, MergesDefaultValues)
{
    builder_.rootObject().addValue<int>("b", 1);
    options_.addOption(gmx::IntegerOption("a").defaultValue(2));
    options_.addOption(gmx::IntegerOption("b").defaultValue(3));
    runTest();
}

TEST_F(TreeValueSupportAdjustTest, OrdersValues)
{
    builder_.rootObject().addValue<int>("a", 1);
    builder_.rootObject().addValue<int>("c", 1);
    builder_.rootObject().addValue<int>("b", 1);
    options_.addOption(gmx::IntegerOption("b").defaultValue(2));
    options_.addOption(gmx::IntegerOption("a").defaultValue(1));
    options_.addOption(gmx::IntegerOption("c").defaultValue(3));
    // TODO: This does not actually test the correct ordering, since the
    // reference data is not currently order-sensitive, but the order can be
    // checked manually from the reference data.
    runTest();
}

} // namespace
