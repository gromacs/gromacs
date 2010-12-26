/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Tests option assignment.
 *
 * In addition to testing gmx::OptionsAssigner, these are the main
 * tests for the classes from basicoptions.h and basicoptionstorage.h (and
 * their base classes) that actually implement the behavior, as well as for the
 * internal implementation of the gmx::Options and gmx::Option classes.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/errorreporting/emptyerrorreporter.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"

namespace
{

TEST(OptionsAssignerTest, HandlesMissingRequiredParameter)
{
    gmx::Options options(NULL, NULL);
    int value = 0;
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").store(&value).required());

    gmx::EmptyErrorReporter errors;
    EXPECT_NE(0, options.finish(&errors));
}

TEST(OptionsAssignerTest, HandlesInvalidMultipleParameter)
{
    gmx::Options options(NULL, NULL);
    std::vector<int> values;
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").storeVector(&values).multiValue());

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.appendValue("1"));
    EXPECT_NE(0, assigner.startOption("p"));
    EXPECT_NE(0, assigner.finish());
}

TEST(OptionsAssignerTest, HandlesMultipleParameter)
{
    gmx::Options options(NULL, NULL);
    std::vector<int> values;
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").storeVector(&values).allowMultiple());

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.appendValue("1"));
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.appendValue("2"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));
    EXPECT_TRUE(options.isSet("p"));
    ASSERT_EQ(2U, values.size());
    EXPECT_EQ(1, values[0]);
    EXPECT_EQ(2, values[1]);
}

TEST(OptionsAssignerTest, HandlesMissingValue)
{
    gmx::Options options(NULL, NULL);
    int value1 = 0, value2 = 0;
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").store(&value1));
    options.addOption(IntegerOption("q").store(&value2));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.startOption("q"));
    EXPECT_EQ(0, assigner.appendValue("2"));
    EXPECT_NE(0, assigner.finish());
}

TEST(OptionsAssignerTest, HandlesExtraValue)
{
    gmx::Options options(NULL, NULL);
    int value1 = 0;
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").store(&value1));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.appendValue("2"));
    EXPECT_NE(0, assigner.appendValue("3"));
    EXPECT_NE(0, assigner.finish());
}

TEST(OptionsAssignerTest, HandlesSubSections)
{
    gmx::Options options(NULL, NULL);
    gmx::Options sub1("section1", NULL);
    gmx::Options sub2("section2", NULL);
    int value = 3;
    int value1 = 1;
    int value2 = 2;
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").store(&value));
    sub1.addOption(IntegerOption("p").store(&value1));
    sub2.addOption(IntegerOption("p").store(&value2));
    options.addSubSection(&sub1);
    options.addSubSection(&sub2);

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.startSubSection("section1"));
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.appendValue("5"));
    EXPECT_EQ(0, assigner.finishSubSection());
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.appendValue("4"));
    EXPECT_EQ(0, assigner.startSubSection("section2"));
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.appendValue("6"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

    EXPECT_EQ(4, value);
    EXPECT_EQ(5, value1);
    EXPECT_EQ(6, value2);
}

TEST(OptionsAssignerTest, HandlesNoStrictSubSections)
{
    gmx::Options options(NULL, NULL);
    gmx::Options sub1("section1", NULL);
    gmx::Options sub2("section2", NULL);
    int pvalue = 3;
    int pvalue1 = 1;
    int qvalue  = 4;
    int pvalue2 = 2;
    int rvalue  = 5;
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").store(&pvalue));
    sub1.addOption(IntegerOption("p").store(&pvalue1));
    sub1.addOption(IntegerOption("q").store(&qvalue));
    sub2.addOption(IntegerOption("p").store(&pvalue2));
    sub2.addOption(IntegerOption("r").store(&rvalue));
    options.addSubSection(&sub1);
    options.addSubSection(&sub2);

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    assigner.setNoStrictSectioning(true);
    EXPECT_EQ(0, assigner.startOption("q"));
    EXPECT_EQ(0, assigner.appendValue("6"));
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.appendValue("7"));
    EXPECT_EQ(0, assigner.startOption("r"));
    EXPECT_EQ(0, assigner.appendValue("8"));
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.appendValue("9"));
    EXPECT_EQ(0, assigner.finishSubSection());
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.appendValue("10"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

    EXPECT_EQ(6, qvalue);
    EXPECT_EQ(7, pvalue1);
    EXPECT_EQ(8, rvalue);
    EXPECT_EQ(9, pvalue2);
    EXPECT_EQ(10, pvalue);
}

TEST(OptionsAssignerTest, HandlesMultipleSources)
{
    gmx::Options options(NULL, NULL);
    int value = -1;
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").store(&value));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.start());
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.appendValue("1"));
    EXPECT_EQ(0, assigner.finish());
    gmx::OptionsAssigner assigner2(&options, &errors);
    EXPECT_EQ(0, assigner2.start());
    EXPECT_EQ(0, assigner2.startOption("p"));
    EXPECT_EQ(0, assigner2.appendValue("2"));
    EXPECT_EQ(0, assigner2.finish());
    EXPECT_EQ(0, options.finish(&errors));

    EXPECT_EQ(2, value);
}


TEST(OptionsAssignerBooleanTest, StoresYesValue)
{
    gmx::Options options(NULL, NULL);
    bool  value = false;
    using gmx::BooleanOption;
    options.addOption(BooleanOption("p").store(&value));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.appendValue("yes"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

    EXPECT_TRUE(value);
}

TEST(OptionsAssignerBooleanTest, SetsBooleanWithoutExplicitValue)
{
    gmx::Options options(NULL, NULL);
    bool value = false;
    using gmx::BooleanOption;
    options.addOption(BooleanOption("p").store(&value));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.startOption("p"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

    EXPECT_TRUE(value);
}

TEST(OptionsAssignerBooleanTest, ClearsBooleanWithPrefixNo)
{
    gmx::Options options(NULL, NULL);
    bool value = true;
    using gmx::BooleanOption;
    options.addOption(BooleanOption("p").store(&value));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    assigner.setAcceptBooleanNoPrefix(true);
    EXPECT_EQ(0, assigner.startOption("nop"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

    EXPECT_FALSE(value);
}

TEST(OptionsAssignerBooleanTest, HandlesBooleanWithPrefixAndValue)
{
    gmx::Options options(NULL, NULL);
    bool value = false;
    using gmx::BooleanOption;
    options.addOption(BooleanOption("p").store(&value));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    assigner.setAcceptBooleanNoPrefix(true);
    EXPECT_EQ(0, assigner.startOption("nop"));
    assigner.appendValue("no");
    int rc = assigner.finish();
    // It's OK to fail, but if it doesn't, it should work.
    if (rc == 0)
    {
        EXPECT_EQ(0, options.finish(&errors));
        EXPECT_TRUE(value);
    }
}


TEST(OptionsAssignerIntegerTest, StoresSingleValue)
{
    gmx::Options options(NULL, NULL);
    int value = 1;
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").store(&value));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    ASSERT_EQ(0, assigner.startOption("p"));
    ASSERT_EQ(0, assigner.appendValue("3"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

    EXPECT_EQ(3, value);
}

TEST(OptionsAssignerIntegerTest, StoresDefaultValue)
{
    gmx::Options options(NULL, NULL);
    int value = -1;
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").store(&value).defaultValue(2));
    EXPECT_EQ(2, value);

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

    EXPECT_EQ(2, value);
}

TEST(OptionsAssignerIntegerTest, StoresToArray)
{
    gmx::Options          options(NULL, NULL);
    int                  *values = NULL;
    int                   nval = -1;
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").storeArray(&values, &nval).multiValue());

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    ASSERT_EQ(0, assigner.startOption("p"));
    ASSERT_EQ(0, assigner.appendValue("-2"));
    ASSERT_EQ(0, assigner.appendValue("1"));
    ASSERT_EQ(0, assigner.appendValue("4"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

    EXPECT_EQ(3, nval);
    EXPECT_EQ(-2, values[0]);
    EXPECT_EQ(1, values[1]);
    EXPECT_EQ(4, values[2]);
    delete[] values;
}

TEST(OptionsAssignerIntegerTest, StoresToVector)
{
    gmx::Options          options(NULL, NULL);
    std::vector<int>      values;
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").storeVector(&values).multiValue());

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    ASSERT_EQ(0, assigner.startOption("p"));
    ASSERT_EQ(0, assigner.appendValue("-2"));
    ASSERT_EQ(0, assigner.appendValue("1"));
    ASSERT_EQ(0, assigner.appendValue("4"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

    EXPECT_EQ(3U, values.size());
    EXPECT_EQ(-2, values[0]);
    EXPECT_EQ(1, values[1]);
    EXPECT_EQ(4, values[2]);
}

TEST(OptionsAssignerIntegerTest, HandlesVectors)
{
    gmx::Options options(NULL, NULL);
    int  vec[3] = {0, 0, 0};
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").store(vec).vector());

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    ASSERT_EQ(0, assigner.startOption("p"));
    ASSERT_EQ(0, assigner.appendValue("-2"));
    ASSERT_EQ(0, assigner.appendValue("1"));
    ASSERT_EQ(0, assigner.appendValue("4"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

    EXPECT_EQ(-2, vec[0]);
    EXPECT_EQ(1, vec[1]);
    EXPECT_EQ(4, vec[2]);
}

TEST(OptionsAssignerIntegerTest, HandlesVectorFromSingleValue)
{
    gmx::Options options(NULL, NULL);
    int  vec[3] = {0, 0, 0};
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").store(vec).vector());

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    ASSERT_EQ(0, assigner.startOption("p"));
    ASSERT_EQ(0, assigner.appendValue("2"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, options.finish(&errors));

    EXPECT_EQ(2, vec[0]);
    EXPECT_EQ(2, vec[1]);
    EXPECT_EQ(2, vec[2]);
}

TEST(OptionsAssignerIntegerTest, HandlesVectorsWithDefaultValue)
{
    gmx::Options options(NULL, NULL);
    int  vec[3] = {3, 2, 1};
    using gmx::IntegerOption;
    options.addOption(IntegerOption("p").store(vec).vector());

    gmx::EmptyErrorReporter errors;
    EXPECT_EQ(0, options.finish(&errors));

    EXPECT_EQ(3, vec[0]);
    EXPECT_EQ(2, vec[1]);
    EXPECT_EQ(1, vec[2]);
}


TEST(OptionsAssignerDoubleTest, StoresSingleValue)
{
    gmx::Options options(NULL, NULL);
    double value = 0.0;
    using gmx::DoubleOption;
    options.addOption(DoubleOption("p").store(&value));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    ASSERT_EQ(0, assigner.startOption("p"));
    ASSERT_EQ(0, assigner.appendValue("2.7"));
    ASSERT_EQ(0, assigner.finish());
    ASSERT_EQ(0, options.finish(&errors));

    EXPECT_DOUBLE_EQ(2.7, value);
}


TEST(OptionsAssignerStringTest, StoresSingleValue)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    using gmx::StringOption;
    options.addOption(StringOption("p").store(&value));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    ASSERT_EQ(0, assigner.startOption("p"));
    ASSERT_EQ(0, assigner.appendValue("value"));
    ASSERT_EQ(0, assigner.finish());
    ASSERT_EQ(0, options.finish(&errors));

    EXPECT_EQ("value", value);
}

TEST(OptionsAssignerStringTest, HandlesEnumValue)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    const char * const     allowed[] = { "none", "test", "value", NULL };
    int                    index = -1;
    using gmx::StringOption;
    options.addOption(StringOption("p").store(&value)
                          .enumValue(allowed)
                          .storeEnumIndex(&index));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    ASSERT_EQ(0, assigner.startOption("p"));
    ASSERT_EQ(0, assigner.appendValue("test"));
    ASSERT_EQ(0, assigner.finish());
    ASSERT_EQ(0, options.finish(&errors));

    EXPECT_EQ("test", value);
    EXPECT_EQ(1, index);
}

TEST(OptionsAssignerStringTest, HandlesIncorrectEnumValue)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    const char * const     allowed[] = { "none", "test", "value", NULL };
    int                    index = -1;
    using gmx::StringOption;
    options.addOption(StringOption("p").store(&value)
                          .enumValue(allowed)
                          .storeEnumIndex(&index));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    ASSERT_EQ(0, assigner.startOption("p"));
    ASSERT_NE(0, assigner.appendValue("unknown"));
}

TEST(OptionsAssignerStringTest, CompletesEnumValue)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    const char * const     allowed[] = { "none", "test", "value", NULL };
    int                    index = -1;
    using gmx::StringOption;
    options.addOption(StringOption("p").store(&value)
                          .enumValue(allowed)
                          .storeEnumIndex(&index));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    ASSERT_EQ(0, assigner.startOption("p"));
    ASSERT_EQ(0, assigner.appendValue("te"));
    ASSERT_EQ(0, assigner.finish());
    ASSERT_EQ(0, options.finish(&errors));

    EXPECT_EQ("test", value);
    EXPECT_EQ(1, index);
}

TEST(OptionsAssignerStringTest, HandlesEnumWithNoValue)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    const char * const     allowed[] = { "none", "test", "value", NULL };
    int                    index = -3;
    using gmx::StringOption;
    options.addOption(StringOption("p").store(&value)
                          .enumValue(allowed)
                          .storeEnumIndex(&index));
    EXPECT_TRUE(value.empty());
    EXPECT_EQ(-1, index);

    gmx::EmptyErrorReporter errors;
    ASSERT_EQ(0, options.finish(&errors));

    EXPECT_TRUE(value.empty());
    EXPECT_EQ(-1, index);
}

TEST(OptionsAssignerStringTest, HandlesEnumDefaultValue)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    const char * const     allowed[] = { "none", "test", "value", NULL };
    int                    index = -1;
    using gmx::StringOption;
    options.addOption(StringOption("p").store(&value)
                          .enumValue(allowed)
                          .defaultValue("test")
                          .storeEnumIndex(&index));
    EXPECT_EQ("test", value);
    EXPECT_EQ(1, index);

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    ASSERT_EQ(0, assigner.finish());
    ASSERT_EQ(0, options.finish(&errors));

    EXPECT_EQ("test", value);
    EXPECT_EQ(1, index);
}

TEST(OptionsAssignerStringTest, HandlesEnumDefaultIndex)
{
    gmx::Options           options(NULL, NULL);
    std::string            value;
    const char * const     allowed[] = { "none", "test", "value", NULL };
    int                    index = -1;
    using gmx::StringOption;
    options.addOption(StringOption("p").store(&value)
                          .enumValue(allowed)
                          .defaultEnumIndex(1)
                          .storeEnumIndex(&index));
    EXPECT_EQ("test", value);
    EXPECT_EQ(1, index);

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&options, &errors);
    ASSERT_EQ(0, assigner.finish());
    ASSERT_EQ(0, options.finish(&errors));

    EXPECT_EQ("test", value);
    EXPECT_EQ(1, index);
}

} // namespace
