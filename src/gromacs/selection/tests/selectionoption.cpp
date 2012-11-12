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
 * Tests handling of selection options.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#include <gtest/gtest.h>

#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionfileoption.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/selection/selectionoptionmanager.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testfilemanager.h"

using gmx::test::TestFileManager;

namespace
{

/********************************************************************
 * Base fixture for tests in this file.
 */

class SelectionOptionTestBase : public ::testing::Test
{
    public:
        SelectionOptionTestBase();

        void setManager();

        gmx::SelectionCollection    sc_;
        gmx::SelectionOptionManager manager_;
        gmx::Options options_;
};

SelectionOptionTestBase::SelectionOptionTestBase()
    : manager_(&sc_), options_(NULL, NULL)
{
    sc_.setReferencePosType("atom");
    sc_.setOutputPosType("atom");
}

void SelectionOptionTestBase::setManager()
{
    setManagerForSelectionOptions(&options_, &manager_);
}


/********************************************************************
 * Tests for SelectionOption
 */

//! Test fixture for gmx::SelectionOption.
typedef SelectionOptionTestBase SelectionOptionTest;

TEST_F(SelectionOptionTest, ParsesSimpleSelection)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(SelectionOption("sel").store(&sel)));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("resname RA RB"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
}


TEST_F(SelectionOptionTest, HandlesDynamicSelectionWhenStaticRequired)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel").store(&sel).onlyStatic()));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_THROW(assigner.appendValue("resname RA RB and x < 5"), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
}


TEST_F(SelectionOptionTest, HandlesTooManySelections)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(SelectionOption("sel").store(&sel)));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("resname RA RB"));
    EXPECT_THROW(assigner.appendValue("resname RB RC"), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
}


TEST_F(SelectionOptionTest, HandlesTooFewSelections)
{
    gmx::Selection sel[2];
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel").store(sel).valueCount(2)));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("resname RA RB"));
    EXPECT_THROW(assigner.finishOption(), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
}


TEST_F(SelectionOptionTest, HandlesAdjuster)
{
    gmx::SelectionList        sel;
    using gmx::SelectionOption;
    gmx::SelectionOptionInfo *info = options_.addOption(
            SelectionOption("sel").storeVector(&sel).multiValue());
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("resname RA RB"));
    EXPECT_NO_THROW(assigner.appendValue("resname RB RC"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
    EXPECT_NO_THROW(info->setValueCount(2));
}


TEST_F(SelectionOptionTest, HandlesDynamicWhenStaticRequiredWithAdjuster)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    gmx::SelectionOptionInfo *info = options_.addOption(
            SelectionOption("sel").store(&sel));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("x < 5"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
    EXPECT_THROW(info->setOnlyStatic(true), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, HandlesTooManySelectionsWithAdjuster)
{
    gmx::SelectionList        sel;
    using gmx::SelectionOption;
    gmx::SelectionOptionInfo *info = options_.addOption(
            SelectionOption("sel").storeVector(&sel).multiValue());
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("resname RA RB"));
    EXPECT_NO_THROW(assigner.appendValue("resname RB RC"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
    EXPECT_THROW(info->setValueCount(1), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, HandlesTooFewSelectionsWithAdjuster)
{
    gmx::SelectionList        sel;
    using gmx::SelectionOption;
    gmx::SelectionOptionInfo *info = options_.addOption(
            SelectionOption("sel").storeVector(&sel).multiValue());
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("resname RA RB"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
    EXPECT_THROW(info->setValueCount(2), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, HandlesDelayedRequiredSelection)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel").store(&sel).required()));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
    ASSERT_NO_THROW(manager_.parseRequestedFromString("resname RA RB"));
    ASSERT_STREQ("resname RA RB", sel.selectionText());
}


TEST_F(SelectionOptionTest, HandlesTooFewDelayedRequiredSelections)
{
    gmx::Selection sel[2];
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel").store(sel).required()
                            .valueCount(2)));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
    EXPECT_THROW(manager_.parseRequestedFromString("resname RA RB"), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, HandlesDelayedOptionalSelection)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(SelectionOption("sel").store(&sel)));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
    ASSERT_NO_THROW(manager_.parseRequestedFromString("resname RA RB"));
    ASSERT_STREQ("resname RA RB", sel.selectionText());
}


TEST_F(SelectionOptionTest, HandlesDelayedSelectionWithAdjuster)
{
    gmx::SelectionList        sel;
    using gmx::SelectionOption;
    gmx::SelectionOptionInfo *info = options_.addOption(
            SelectionOption("sel").storeVector(&sel).valueCount(3));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
    EXPECT_NO_THROW(info->setValueCount(2));
    EXPECT_NO_THROW(manager_.parseRequestedFromString("resname RA RB; resname RB RC"));
}


/********************************************************************
 * Tests for SelectionFileOption
 */

class SelectionFileOptionTest : public SelectionOptionTestBase
{
    public:
        SelectionFileOptionTest();
};

SelectionFileOptionTest::SelectionFileOptionTest()
{
    options_.addOption(gmx::SelectionFileOption("sf"));
}


TEST_F(SelectionFileOptionTest, HandlesSingleSelectionOptionFromFile)
{
    gmx::SelectionList sel;
    gmx::SelectionList reqsel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel").storeVector(&sel).multiValue()));
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("reqsel").storeVector(&reqsel)
                            .multiValue().required()));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startOption("sf"));
    EXPECT_NO_THROW(assigner.appendValue(TestFileManager::getInputFilePath("selfile.dat")));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());

    // These should match the contents of selfile.dat
    ASSERT_EQ(2U, sel.size());
    EXPECT_STREQ("resname RA RB", sel[0].selectionText());
    EXPECT_STREQ("resname RB RC", sel[1].selectionText());
    ASSERT_EQ(0U, reqsel.size());
}


TEST_F(SelectionFileOptionTest, HandlesTwoSeparateSelectionOptions)
{
    gmx::SelectionList sel1;
    gmx::SelectionList sel2;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel1").storeVector(&sel1).multiValue()));
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel2").storeVector(&sel2).multiValue()));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    std::string          value(TestFileManager::getInputFilePath("selfile.dat"));
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel1"));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startOption("sf"));
    EXPECT_NO_THROW(assigner.appendValue(value));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startOption("sel2"));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startOption("sf"));
    EXPECT_NO_THROW(assigner.appendValue(value));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());

    // These should match the contents of selfile.dat
    ASSERT_EQ(2U, sel1.size());
    EXPECT_STREQ("resname RA RB", sel1[0].selectionText());
    EXPECT_STREQ("resname RB RC", sel1[1].selectionText());
    ASSERT_EQ(2U, sel2.size());
    EXPECT_STREQ("resname RA RB", sel2[0].selectionText());
    EXPECT_STREQ("resname RB RC", sel2[1].selectionText());
}


TEST_F(SelectionFileOptionTest, HandlesTwoSelectionOptionsFromSingleFile)
{
    gmx::SelectionList sel1;
    gmx::SelectionList sel2;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel1").storeVector(&sel1)));
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel2").storeVector(&sel2)));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    std::string          value(TestFileManager::getInputFilePath("selfile.dat"));
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel1"));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startOption("sel2"));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startOption("sf"));
    EXPECT_NO_THROW(assigner.appendValue(value));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());

    // These should match the contents of selfile.dat
    ASSERT_EQ(1U, sel1.size());
    EXPECT_STREQ("resname RA RB", sel1[0].selectionText());
    ASSERT_EQ(1U, sel2.size());
    EXPECT_STREQ("resname RB RC", sel2[0].selectionText());
}


TEST_F(SelectionFileOptionTest, HandlesRequiredOptionFromFile)
{
    gmx::SelectionList sel;
    gmx::SelectionList optsel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel").storeVector(&sel)
                            .multiValue().required()));
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("optsel").storeVector(&optsel)
                            .multiValue()));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sf"));
    EXPECT_NO_THROW(assigner.appendValue(TestFileManager::getInputFilePath("selfile.dat")));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.startOption("optsel"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
    EXPECT_NO_THROW(manager_.parseRequestedFromString("resname RC RD"));

    // These should match the contents of selfile.dat
    ASSERT_EQ(2U, sel.size());
    EXPECT_STREQ("resname RA RB", sel[0].selectionText());
    EXPECT_STREQ("resname RB RC", sel[1].selectionText());
    ASSERT_EQ(1U, optsel.size());
    EXPECT_STREQ("resname RC RD", optsel[0].selectionText());
}


// TODO: Is this the best possible behavior, or should it error out?
TEST_F(SelectionFileOptionTest, HandlesRequiredOptionFromFileWithOtherOptionSet)
{
    gmx::SelectionList sel1;
    gmx::SelectionList sel2;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel1").storeVector(&sel1)
                            .multiValue().required()));
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel2").storeVector(&sel2)
                            .multiValue().required()));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.startOption("sel1"));
    EXPECT_NO_THROW(assigner.appendValue("resname RC RD"));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startOption("sf"));
    EXPECT_NO_THROW(assigner.appendValue(TestFileManager::getInputFilePath("selfile.dat")));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());

    // These should match the contents of selfile.dat
    ASSERT_EQ(2U, sel2.size());
    EXPECT_STREQ("resname RA RB", sel2[0].selectionText());
    EXPECT_STREQ("resname RB RC", sel2[1].selectionText());
    ASSERT_EQ(1U, sel1.size());
    EXPECT_STREQ("resname RC RD", sel1[0].selectionText());
}


TEST_F(SelectionFileOptionTest, HandlesTwoRequiredOptionsFromSingleFile)
{
    gmx::SelectionList sel1;
    gmx::SelectionList sel2;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel1").storeVector(&sel1).required()));
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel2").storeVector(&sel2).required()));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    std::string          value(TestFileManager::getInputFilePath("selfile.dat"));
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sf"));
    EXPECT_NO_THROW(assigner.appendValue(value));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());

    // These should match the contents of selfile.dat
    ASSERT_EQ(1U, sel1.size());
    EXPECT_STREQ("resname RA RB", sel1[0].selectionText());
    ASSERT_EQ(1U, sel2.size());
    EXPECT_STREQ("resname RB RC", sel2[0].selectionText());
}


TEST_F(SelectionFileOptionTest, GivesErrorWithNoFile)
{
    gmx::SelectionList sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel").storeVector(&sel).multiValue()));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startOption("sf"));
    EXPECT_THROW(assigner.finishOption(), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
}


TEST_F(SelectionFileOptionTest, GivesErrorWithNonExistentFile)
{
    gmx::SelectionList sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel").storeVector(&sel).multiValue()));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startOption("sf"));
    // TODO: Should this be changed to an InvalidInputError?
    EXPECT_THROW(assigner.appendValue("nonexistentfile"), gmx::FileIOError);
    EXPECT_THROW(assigner.appendValue(TestFileManager::getInputFilePath("selfile.dat")),
                 gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
}


TEST_F(SelectionFileOptionTest, GivesErrorWithMultipleFiles)
{
    gmx::SelectionList sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(options_.addOption(
                        SelectionOption("sel").storeVector(&sel).multiValue()));
    setManager();

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.finishOption());
    ASSERT_NO_THROW(assigner.startOption("sf"));
    EXPECT_NO_THROW(assigner.appendValue(TestFileManager::getInputFilePath("selfile.dat")));
    EXPECT_THROW(assigner.appendValue("nonexistentfile"), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(options_.finish());
}

} // namespace
