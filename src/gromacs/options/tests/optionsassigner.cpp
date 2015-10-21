/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Tests option assignment.
 *
 * In addition to testing gmx::OptionsAssigner, these are the main
 * tests for the classes from basicoptions.h and basicoptionstorage.h (and
 * their base classes) that actually implement the behavior, as well as for the
 * internal implementation of the gmx::Options and gmx::AbstractOptionStorage
 * classes.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "gromacs/options/optionsassigner.h"

#include <limits>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"

namespace
{

/********************************************************************
 * General assignment tests
 */

TEST(OptionsAssignerTest, HandlesMissingRequiredParameter)
{
    gmx::Options options(NULL, NULL);
    int          value = 0;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(
                            IntegerOption("p").store(&value).required()));

    EXPECT_THROW(options.finish(), gmx::InvalidInputError);
}

TEST(OptionsAssignerTest, HandlesRequiredParameterWithDefaultValue)
{
    gmx::Options options(NULL, NULL);
    int          value = 0;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(
                            IntegerOption("p").store(&value).required()
                                .defaultValue(1)));

    EXPECT_NO_THROW(options.finish());
    EXPECT_EQ(1, value);
}

TEST(OptionsAssignerTest, HandlesInvalidMultipleParameter)
{
    gmx::Options     options(NULL, NULL);
    std::vector<int> values;
    bool             bIsSet = false;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(
                            IntegerOption("p")
                                .storeVector(&values).storeIsSet(&bIsSet)
                                .multiValue()));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("1"));
    ASSERT_NO_THROW(assigner.finishOption());
    EXPECT_THROW(assigner.startOption("p"), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_TRUE(bIsSet);
    ASSERT_EQ(1U, values.size());
    EXPECT_EQ(1, values[0]);
}

TEST(OptionsAssignerTest, HandlesMultipleParameter)
{
    gmx::Options     options(NULL, NULL);
    std::vector<int> values;
    bool             bIsSet = false;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(
                            IntegerOption("p")
                                .storeVector(&values).storeIsSet(&bIsSet)
                                .allowMultiple()));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("p"));
    ASSERT_NO_THROW_GMX(assigner.appendValue("1"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    ASSERT_NO_THROW_GMX(assigner.startOption("p"));
    ASSERT_NO_THROW_GMX(assigner.appendValue("2"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options.finish());

    EXPECT_TRUE(bIsSet);
    ASSERT_EQ(2U, values.size());
    EXPECT_EQ(1, values[0]);
    EXPECT_EQ(2, values[1]);
}

TEST(OptionsAssignerTest, HandlesMissingValue)
{
    gmx::Options options(NULL, NULL);
    int          value1 = 0, value2 = 0;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(IntegerOption("p").store(&value1)));
    ASSERT_NO_THROW(options.addOption(IntegerOption("q").store(&value2)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    EXPECT_THROW(assigner.finishOption(), gmx::InvalidInputError);
    ASSERT_NO_THROW(assigner.startOption("q"));
    ASSERT_NO_THROW(assigner.appendValue("2"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(2, value2);
}

TEST(OptionsAssignerTest, HandlesExtraValue)
{
    gmx::Options options(NULL, NULL);
    int          value1 = 0;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(IntegerOption("p").store(&value1)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("2"));
    EXPECT_THROW(assigner.appendValue("3"), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(0, value1);
}

TEST(OptionsAssignerTest, HandlesGroups)
{
    gmx::Options            options(NULL, NULL);
    gmx::IOptionsContainer &group1 = options.addGroup();
    gmx::IOptionsContainer &group2 = options.addGroup();
    int                     value  = 3;
    int                     value1 = 1;
    int                     value2 = 2;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(IntegerOption("p").store(&value)));
    ASSERT_NO_THROW(group1.addOption(IntegerOption("q").store(&value1)));
    ASSERT_NO_THROW(group2.addOption(IntegerOption("r").store(&value2)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    EXPECT_NO_THROW(assigner.appendValue("5"));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startOption("q"));
    EXPECT_NO_THROW(assigner.appendValue("4"));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startOption("r"));
    EXPECT_NO_THROW(assigner.appendValue("6"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(5, value);
    EXPECT_EQ(4, value1);
    EXPECT_EQ(6, value2);
}

TEST(OptionsAssignerTest, HandlesSubSections)
{
    gmx::Options options(NULL, NULL);
    gmx::Options sub1("section1", NULL);
    gmx::Options sub2("section2", NULL);
    int          value  = 3;
    int          value1 = 1;
    int          value2 = 2;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addSubSection(&sub1));
    ASSERT_NO_THROW(options.addSubSection(&sub2));
    ASSERT_NO_THROW(options.addOption(IntegerOption("p").store(&value)));
    ASSERT_NO_THROW(sub1.addOption(IntegerOption("p").store(&value1)));
    ASSERT_NO_THROW(sub2.addOption(IntegerOption("p").store(&value2)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startSubSection("section1"));
    ASSERT_NO_THROW(assigner.startOption("p"));
    EXPECT_NO_THROW(assigner.appendValue("5"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finishSubSection());
    ASSERT_NO_THROW(assigner.startOption("p"));
    EXPECT_NO_THROW(assigner.appendValue("4"));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startSubSection("section2"));
    ASSERT_NO_THROW(assigner.startOption("p"));
    EXPECT_NO_THROW(assigner.appendValue("6"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finishSubSection());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(4, value);
    EXPECT_EQ(5, value1);
    EXPECT_EQ(6, value2);
}

TEST(OptionsAssignerTest, HandlesMultipleSources)
{
    gmx::Options options(NULL, NULL);
    int          value = -1;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(IntegerOption("p").store(&value)));

    {
        gmx::OptionsAssigner assigner(&options);
        EXPECT_NO_THROW(assigner.start());
        ASSERT_NO_THROW(assigner.startOption("p"));
        EXPECT_NO_THROW(assigner.appendValue("1"));
        EXPECT_NO_THROW(assigner.finishOption());
        EXPECT_NO_THROW(assigner.finish());
    }
    {
        gmx::OptionsAssigner assigner2(&options);
        EXPECT_NO_THROW(assigner2.start());
        ASSERT_NO_THROW(assigner2.startOption("p"));
        EXPECT_NO_THROW(assigner2.appendValue("2"));
        EXPECT_NO_THROW(assigner2.finishOption());
        EXPECT_NO_THROW(assigner2.finish());
    }
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(2, value);
}


/********************************************************************
 * Tests for boolean assignment
 */

TEST(OptionsAssignerBooleanTest, StoresYesValue)
{
    gmx::Options options(NULL, NULL);
    bool         value = false;
    using gmx::BooleanOption;
    ASSERT_NO_THROW(options.addOption(BooleanOption("p").store(&value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    EXPECT_NO_THROW(assigner.appendValue("yes"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_TRUE(value);
}

TEST(OptionsAssignerBooleanTest, SetsBooleanWithoutExplicitValue)
{
    gmx::Options options(NULL, NULL);
    bool         value = false;
    using gmx::BooleanOption;
    ASSERT_NO_THROW(options.addOption(BooleanOption("p").store(&value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_TRUE(value);
}

TEST(OptionsAssignerBooleanTest, ClearsBooleanWithPrefixNo)
{
    gmx::Options options(NULL, NULL);
    bool         value = true;
    using gmx::BooleanOption;
    ASSERT_NO_THROW(options.addOption(BooleanOption("p").store(&value)));

    gmx::OptionsAssigner assigner(&options);
    assigner.setAcceptBooleanNoPrefix(true);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("nop"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_FALSE(value);
}

TEST(OptionsAssignerBooleanTest, HandlesBooleanWithPrefixAndValue)
{
    gmx::Options options(NULL, NULL);
    bool         value = false;
    using gmx::BooleanOption;
    ASSERT_NO_THROW(options.addOption(BooleanOption("p").store(&value)));

    gmx::OptionsAssigner assigner(&options);
    assigner.setAcceptBooleanNoPrefix(true);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("nop"));
    // It's OK to fail, but if it doesn't, it should work.
    try
    {
        assigner.appendValue("no");
        assigner.finishOption();
        EXPECT_NO_THROW(assigner.finish());
        EXPECT_TRUE(value);
    }
    catch (const gmx::InvalidInputError &)
    {
    }
}


/********************************************************************
 * Tests for integer assignment
 *
 * These tests also contain tests for general default value handling.
 */

TEST(OptionsAssignerIntegerTest, StoresSingleValue)
{
    gmx::Options options(NULL, NULL);
    int          value = 1;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(IntegerOption("p").store(&value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("3"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(3, value);
}

TEST(OptionsAssignerIntegerTest, HandlesEmptyValue)
{
    gmx::Options options(NULL, NULL);
    int          value = 1;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(IntegerOption("p").store(&value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    EXPECT_THROW(assigner.appendValue(""), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(1, value);
}

TEST(OptionsAssignerIntegerTest, HandlesInvalidValue)
{
    gmx::Options options(NULL, NULL);
    int          value = 1;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(IntegerOption("p").store(&value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    EXPECT_THROW(assigner.appendValue("2abc"), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(1, value);
}

TEST(OptionsAssignerIntegerTest, HandlesOverflow)
{
    gmx::Options options(NULL, NULL);
    int          value = 1;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(IntegerOption("p").store(&value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    std::string overflowValue(
            gmx::formatString("%d0000", std::numeric_limits<int>::max()));
    EXPECT_THROW(assigner.appendValue(overflowValue), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(1, value);
}

TEST(OptionsAssignerIntegerTest, StoresDefaultValue)
{
    gmx::Options options(NULL, NULL);
    int          value = -1;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(
                            IntegerOption("p").store(&value).defaultValue(2)));
    EXPECT_EQ(2, value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(2, value);
}

TEST(OptionsAssignerIntegerTest, StoresDefaultValueIfSet)
{
    gmx::Options options(NULL, NULL);
    int          value = -1;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(
                            IntegerOption("p").store(&value).defaultValueIfSet(2)));
    EXPECT_EQ(-1, value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(2, value);
}

TEST(OptionsAssignerIntegerTest, HandlesDefaultValueIfSetWhenNotSet)
{
    gmx::Options options(NULL, NULL);
    int          value = -1;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(
                            IntegerOption("p").store(&value).defaultValueIfSet(2)));
    EXPECT_EQ(-1, value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(-1, value);
}

TEST(OptionsAssignerIntegerTest, HandlesBothDefaultValues)
{
    gmx::Options options(NULL, NULL);
    int          value = -1;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(
                            IntegerOption("p").store(&value)
                                .defaultValue(1).defaultValueIfSet(2)));
    EXPECT_EQ(1, value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(2, value);
}

TEST(OptionsAssignerIntegerTest, StoresToVector)
{
    gmx::Options          options(NULL, NULL);
    std::vector<int>      values;
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(
                            IntegerOption("p").storeVector(&values).multiValue()));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("-2"));
    ASSERT_NO_THROW(assigner.appendValue("1"));
    ASSERT_NO_THROW(assigner.appendValue("4"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(3U, values.size());
    EXPECT_EQ(-2, values[0]);
    EXPECT_EQ(1, values[1]);
    EXPECT_EQ(4, values[2]);
}

TEST(OptionsAssignerIntegerTest, HandlesVectors)
{
    gmx::Options options(NULL, NULL);
    int          vec[3] = {0, 0, 0};
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(IntegerOption("p").store(vec).vector()));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("-2"));
    ASSERT_NO_THROW(assigner.appendValue("1"));
    ASSERT_NO_THROW(assigner.appendValue("4"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(-2, vec[0]);
    EXPECT_EQ(1, vec[1]);
    EXPECT_EQ(4, vec[2]);
}

TEST(OptionsAssignerIntegerTest, HandlesVectorFromSingleValue)
{
    gmx::Options options(NULL, NULL);
    int          vec[3] = {0, 0, 0};
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(IntegerOption("p").store(vec).vector()));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("2"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(2, vec[0]);
    EXPECT_EQ(2, vec[1]);
    EXPECT_EQ(2, vec[2]);
}

TEST(OptionsAssignerIntegerTest, HandlesVectorsWithDefaultValue)
{
    gmx::Options options(NULL, NULL);
    int          vec[3] = {3, 2, 1};
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(IntegerOption("p").store(vec).vector()));

    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(3, vec[0]);
    EXPECT_EQ(2, vec[1]);
    EXPECT_EQ(1, vec[2]);
}

TEST(OptionsAssignerIntegerTest, HandlesVectorsWithDefaultValueWithInvalidAssignment)
{
    gmx::Options     options(NULL, NULL);
    int              vec[3] = {3, 2, 1};
    std::vector<int> vec2(vec, vec+3);
    using gmx::IntegerOption;
    ASSERT_NO_THROW(options.addOption(IntegerOption("p").store(vec)
                                          .storeVector(&vec2).vector()));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    EXPECT_NO_THROW(assigner.appendValue("1"));
    EXPECT_NO_THROW(assigner.appendValue("3"));
    EXPECT_THROW(assigner.finishOption(), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(3, vec[0]);
    EXPECT_EQ(2, vec[1]);
    EXPECT_EQ(1, vec[2]);
    ASSERT_EQ(3U, vec2.size());
    EXPECT_EQ(3, vec2[0]);
    EXPECT_EQ(2, vec2[1]);
    EXPECT_EQ(1, vec2[2]);
}


/********************************************************************
 * Tests for double assignment
 */

TEST(OptionsAssignerDoubleTest, StoresSingleValue)
{
    gmx::Options options(NULL, NULL);
    double       value = 0.0;
    using gmx::DoubleOption;
    ASSERT_NO_THROW(options.addOption(DoubleOption("p").store(&value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("2.7"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_DOUBLE_EQ(2.7, value);
}

TEST(OptionsAssignerDoubleTest, HandlesEmptyValue)
{
    gmx::Options options(NULL, NULL);
    double       value = 1.0;
    using gmx::DoubleOption;
    ASSERT_NO_THROW(options.addOption(DoubleOption("p").store(&value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    EXPECT_THROW(assigner.appendValue(""), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_DOUBLE_EQ(1.0, value);
}


/********************************************************************
 * Tests for string assignment
 */

//! Set of allowed values for enum option tests.
const char *const c_allowed[] = { "none", "test", "value" };

TEST(OptionsAssignerStringTest, StoresSingleValue)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    using gmx::StringOption;
    ASSERT_NO_THROW(options.addOption(StringOption("p").store(&value)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("value"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ("value", value);
}

TEST(OptionsAssignerStringTest, HandlesEnumValue)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    using gmx::StringOption;
    ASSERT_NO_THROW(options.addOption(
                            StringOption("p").store(&value)
                                .enumValue(c_allowed)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("test"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ("test", value);
}

TEST(OptionsAssignerStringTest, HandlesEnumValueFromNullTerminatedArray)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    const char * const     allowed[] = { "none", "test", "value", NULL };
    using gmx::StringOption;
    ASSERT_NO_THROW(options.addOption(
                            StringOption("p").store(&value)
                                .enumValueFromNullTerminatedArray(allowed)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("value"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ("value", value);
}

TEST(OptionsAssignerStringTest, HandlesIncorrectEnumValue)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    using gmx::StringOption;
    ASSERT_NO_THROW(options.addOption(
                            StringOption("p").store(&value)
                                .enumValue(c_allowed)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_THROW(assigner.appendValue("unknown"), gmx::InvalidInputError);
}

TEST(OptionsAssignerStringTest, CompletesEnumValue)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    using gmx::StringOption;
    ASSERT_NO_THROW(options.addOption(
                            StringOption("p").store(&value)
                                .enumValue(c_allowed)));

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("te"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ("test", value);
}

TEST(OptionsAssignerStringTest, HandlesEnumWithNoValue)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    using gmx::StringOption;
    ASSERT_NO_THROW(options.addOption(
                            StringOption("p").store(&value)
                                .enumValue(c_allowed)));
    EXPECT_TRUE(value.empty());

    ASSERT_NO_THROW(options.finish());

    EXPECT_TRUE(value.empty());
}

TEST(OptionsAssignerStringTest, HandlesEnumDefaultValue)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    using gmx::StringOption;
    ASSERT_NO_THROW(options.addOption(
                            StringOption("p").store(&value)
                                .enumValue(c_allowed).defaultValue("test")));
    EXPECT_EQ("test", value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ("test", value);
}

TEST(OptionsAssignerStringTest, HandlesEnumDefaultValueFromVariable)
{
    gmx::Options           options(NULL, NULL);
    std::string            value("test");
    using gmx::StringOption;
    ASSERT_NO_THROW(options.addOption(
                            StringOption("p").store(&value)
                                .enumValue(c_allowed)));
    EXPECT_EQ("test", value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ("test", value);
}

TEST(OptionsAssignerStringTest, HandlesEnumDefaultValueFromVector)
{
    gmx::Options             options(NULL, NULL);
    std::vector<std::string> value;
    value.push_back("test");
    value.push_back("value");
    using gmx::StringOption;
    ASSERT_NO_THROW(options.addOption(
                            StringOption("p").storeVector(&value).valueCount(2)
                                .enumValue(c_allowed)));
    ASSERT_EQ(2U, value.size());
    EXPECT_EQ("test", value[0]);
    EXPECT_EQ("value", value[1]);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    ASSERT_EQ(2U, value.size());
    EXPECT_EQ("test", value[0]);
    EXPECT_EQ("value", value[1]);
}


/********************************************************************
 * Tests for enum options
 */

//! Enum type for EnumOption tests.
enum TestEnum
{
    etestNone,
    etestTest,
    etestValue,
    etestNR
};

TEST(OptionsAssignerEnumTest, StoresSingleValue)
{
    gmx::Options           options(NULL, NULL);
    TestEnum               value     = etestNone;
    using gmx::EnumOption;
    ASSERT_NO_THROW(options.addOption(
                            EnumOption<TestEnum>("p").store(&value)
                                .enumValue(c_allowed)));
    EXPECT_EQ(etestNone, value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("test"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(etestTest, value);
}

TEST(OptionsAssignerEnumTest, StoresVectorValues)
{
    gmx::Options           options(NULL, NULL);
    std::vector<TestEnum>  values;
    using gmx::EnumOption;
    ASSERT_NO_THROW(options.addOption(
                            EnumOption<TestEnum>("p").storeVector(&values)
                                .multiValue().enumValue(c_allowed)));
    EXPECT_TRUE(values.empty());

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("p"));
    ASSERT_NO_THROW(assigner.appendValue("test"));
    ASSERT_NO_THROW(assigner.appendValue("value"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    ASSERT_EQ(2U, values.size());
    EXPECT_EQ(etestTest, values[0]);
    EXPECT_EQ(etestValue, values[1]);
}

TEST(OptionsAssignerEnumTest, HandlesInitialValueOutOfRange)
{
    gmx::Options           options(NULL, NULL);
    TestEnum               value     = etestNR;
    using gmx::EnumOption;
    ASSERT_NO_THROW(options.addOption(
                            EnumOption<TestEnum>("p").store(&value)
                                .enumValue(c_allowed)));
    EXPECT_EQ(etestNR, value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(etestNR, value);
}

TEST(OptionsAssignerEnumTest, HandlesEnumDefaultValue)
{
    gmx::Options           options(NULL, NULL);
    TestEnum               value     = etestNone;
    using gmx::EnumOption;
    ASSERT_NO_THROW(options.addOption(
                            EnumOption<TestEnum>("p").store(&value)
                                .enumValue(c_allowed).defaultValue(etestTest)));
    EXPECT_EQ(etestTest, value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(etestTest, value);
}

TEST(OptionsAssignerEnumTest, HandlesEnumDefaultValueFromVariable)
{
    gmx::Options           options(NULL, NULL);
    TestEnum               value     = etestTest;
    using gmx::EnumOption;
    ASSERT_NO_THROW(options.addOption(
                            EnumOption<TestEnum>("p").store(&value)
                                .enumValue(c_allowed)));
    EXPECT_EQ(etestTest, value);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    EXPECT_EQ(etestTest, value);
}

TEST(OptionsAssignerEnumTest, HandlesEnumDefaultValueFromVector)
{
    gmx::Options             options(NULL, NULL);
    std::vector<TestEnum>    value;
    value.push_back(etestNone);
    value.push_back(etestTest);
    using gmx::EnumOption;
    ASSERT_NO_THROW(options.addOption(
                            EnumOption<TestEnum>("p").storeVector(&value).valueCount(2)
                                .enumValue(c_allowed)));
    ASSERT_EQ(2U, value.size());
    EXPECT_EQ(etestNone, value[0]);
    EXPECT_EQ(etestTest, value[1]);

    gmx::OptionsAssigner assigner(&options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options.finish());

    ASSERT_EQ(2U, value.size());
    EXPECT_EQ(etestNone, value[0]);
    EXPECT_EQ(etestTest, value[1]);
}

} // namespace
