/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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
 * Tests handling of selection options.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/selection/selectionoption.h"

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionfileoption.h"
#include "gromacs/selection/selectionoptionmanager.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "toputils.h"

using gmx::test::TestFileManager;

namespace gmx
{
namespace test
{
namespace
{

/********************************************************************
 * Base fixture for tests in this file.
 */

class SelectionOptionTestBase : public ::testing::Test
{
public:
    SelectionOptionTestBase();

    void loadTopology(const char* filename);

    gmx::SelectionCollection    sc_;
    gmx::SelectionOptionManager manager_;
    gmx::Options                options_;

private:
    gmx::test::TopologyManager topManager_;
};

SelectionOptionTestBase::SelectionOptionTestBase() : manager_(&sc_)
{
    options_.addManager(&manager_);
    sc_.setReferencePosType("atom");
    sc_.setOutputPosType("atom");
}

void SelectionOptionTestBase::loadTopology(const char* filename)
{
    topManager_.loadTopology(filename);

    ASSERT_NO_THROW_GMX(sc_.setTopology(topManager_.topology(), -1));
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
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").store(&sel)));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("resname RA RB"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    ASSERT_TRUE(sel.isValid());
}


TEST_F(SelectionOptionTest, HandlesDynamicSelectionWhenStaticRequired)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").store(&sel).onlyStatic()));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_THROW_GMX(assigner.appendValue("resname RA RB and x < 5"), gmx::InvalidInputError);
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
}


TEST_F(SelectionOptionTest, HandlesNonAtomicSelectionWhenAtomsRequired)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").store(&sel).onlyAtoms()));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("res_cog of resname RA RB"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW_GMX(sc_.compile(), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, ChecksForSortedAtomsWhenRequired)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").store(&sel).onlySortedAtoms()));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("resnr 1 2 permute 2 1"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW_GMX(sc_.compile(), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, ChecksEmptySelections)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").store(&sel)));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("none"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

    EXPECT_THROW_GMX(sc_.compile(), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, ChecksEmptyDelayedSelections)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").store(&sel)));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
    ASSERT_NO_THROW_GMX(manager_.parseRequestedFromString("none"));

    EXPECT_THROW_GMX(sc_.compile(), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, HandlesTooManySelections)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").store(&sel)));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("resname RA RB"));
    EXPECT_THROW_GMX(assigner.appendValue("resname RB RC"), gmx::InvalidInputError);
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
}


TEST_F(SelectionOptionTest, HandlesTooFewSelections)
{
    gmx::Selection sel[2];
    using gmx::SelectionOption;
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").store(sel).valueCount(2)));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("resname RA RB"));
    EXPECT_THROW_GMX(assigner.finishOption(), gmx::InvalidInputError);
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
}


TEST_F(SelectionOptionTest, HandlesDefaultSelectionText)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    options_.addOption(SelectionOption("sel").store(&sel).defaultSelectionText("all"));

    EXPECT_NO_THROW_GMX(options_.finish());

    ASSERT_TRUE(sel.isValid());

    EXPECT_NO_THROW_GMX(sc_.setTopology(nullptr, 10));
    EXPECT_NO_THROW_GMX(sc_.compile());

    EXPECT_STREQ("all", sel.selectionText());
}


TEST_F(SelectionOptionTest, HandlesAdjuster)
{
    gmx::SelectionList sel;
    using gmx::SelectionOption;
    gmx::SelectionOptionInfo* info =
            options_.addOption(SelectionOption("sel").storeVector(&sel).multiValue());

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("resname RA RB"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("resname RB RC"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
    EXPECT_NO_THROW_GMX(info->setValueCount(2));
}


TEST_F(SelectionOptionTest, HandlesDynamicWhenStaticRequiredWithAdjuster)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    gmx::SelectionOptionInfo* info = options_.addOption(SelectionOption("sel").store(&sel));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("x < 5"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
    EXPECT_THROW_GMX(info->setOnlyStatic(true), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, HandlesTooManySelectionsWithAdjuster)
{
    gmx::SelectionList sel;
    using gmx::SelectionOption;
    gmx::SelectionOptionInfo* info =
            options_.addOption(SelectionOption("sel").storeVector(&sel).multiValue());

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("resname RA RB"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("resname RB RC"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
    EXPECT_THROW_GMX(info->setValueCount(1), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, HandlesTooFewSelectionsWithAdjuster)
{
    gmx::SelectionList sel;
    using gmx::SelectionOption;
    gmx::SelectionOptionInfo* info =
            options_.addOption(SelectionOption("sel").storeVector(&sel).multiValue());

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("resname RA RB"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
    EXPECT_THROW_GMX(info->setValueCount(2), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, HandlesDelayedRequiredSelection)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").store(&sel).required()));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
    ASSERT_NO_THROW_GMX(manager_.parseRequestedFromString("resname RA RB"));
    ASSERT_STREQ("resname RA RB", sel.selectionText());
}


TEST_F(SelectionOptionTest, HandlesTooFewDelayedRequiredSelections)
{
    gmx::Selection sel[2];
    using gmx::SelectionOption;
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").store(sel).required().valueCount(2)));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
    EXPECT_THROW_GMX(manager_.parseRequestedFromString("resname RA RB"), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, HandlesDelayedOptionalSelection)
{
    gmx::Selection sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").store(&sel)));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
    ASSERT_NO_THROW_GMX(manager_.parseRequestedFromString("resname RA RB"));
    ASSERT_STREQ("resname RA RB", sel.selectionText());
}


TEST_F(SelectionOptionTest, HandlesDelayedSelectionWithAdjuster)
{
    gmx::SelectionList sel;
    using gmx::SelectionOption;
    gmx::SelectionOptionInfo* info =
            options_.addOption(SelectionOption("sel").storeVector(&sel).valueCount(3));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
    EXPECT_NO_THROW_GMX(info->setValueCount(2));
    EXPECT_NO_THROW_GMX(manager_.parseRequestedFromString("resname RA RB; resname RB RC"));
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
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").storeVector(&sel).multiValue()));
    ASSERT_NO_THROW_GMX(
            options_.addOption(SelectionOption("reqsel").storeVector(&reqsel).multiValue().required()));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    ASSERT_NO_THROW_GMX(assigner.startOption("sf"));
    EXPECT_NO_THROW_GMX(assigner.appendValue(TestFileManager::getInputFilePath("selfile.dat").string()));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

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
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel1").storeVector(&sel1).multiValue()));
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel2").storeVector(&sel2).multiValue()));

    gmx::OptionsAssigner assigner(&options_);
    std::string          value(TestFileManager::getInputFilePath("selfile.dat").string());
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel1"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    ASSERT_NO_THROW_GMX(assigner.startOption("sf"));
    EXPECT_NO_THROW_GMX(assigner.appendValue(value));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel2"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    ASSERT_NO_THROW_GMX(assigner.startOption("sf"));
    EXPECT_NO_THROW_GMX(assigner.appendValue(value));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

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
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel1").storeVector(&sel1)));
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel2").storeVector(&sel2)));

    gmx::OptionsAssigner assigner(&options_);
    std::string          value(TestFileManager::getInputFilePath("selfile.dat").string());
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel1"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel2"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    ASSERT_NO_THROW_GMX(assigner.startOption("sf"));
    EXPECT_NO_THROW_GMX(assigner.appendValue(value));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

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
    ASSERT_NO_THROW_GMX(
            options_.addOption(SelectionOption("sel").storeVector(&sel).multiValue().required()));
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("optsel").storeVector(&optsel).multiValue()));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sf"));
    EXPECT_NO_THROW_GMX(assigner.appendValue(TestFileManager::getInputFilePath("selfile.dat").string()));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.startOption("optsel"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
    EXPECT_NO_THROW_GMX(manager_.parseRequestedFromString("resname RC RD"));

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
    ASSERT_NO_THROW_GMX(
            options_.addOption(SelectionOption("sel1").storeVector(&sel1).multiValue().required()));
    ASSERT_NO_THROW_GMX(
            options_.addOption(SelectionOption("sel2").storeVector(&sel2).multiValue().required()));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    EXPECT_NO_THROW_GMX(assigner.startOption("sel1"));
    EXPECT_NO_THROW_GMX(assigner.appendValue("resname RC RD"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    ASSERT_NO_THROW_GMX(assigner.startOption("sf"));
    EXPECT_NO_THROW_GMX(assigner.appendValue(TestFileManager::getInputFilePath("selfile.dat").string()));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

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
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel1").storeVector(&sel1).required()));
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel2").storeVector(&sel2).required()));

    gmx::OptionsAssigner assigner(&options_);
    std::string          value(TestFileManager::getInputFilePath("selfile.dat").string());
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sf"));
    EXPECT_NO_THROW_GMX(assigner.appendValue(value));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());

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
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").storeVector(&sel).multiValue()));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    ASSERT_NO_THROW_GMX(assigner.startOption("sf"));
    EXPECT_THROW_GMX(assigner.finishOption(), gmx::InvalidInputError);
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
}


TEST_F(SelectionFileOptionTest, GivesErrorWithNonExistentFile)
{
    gmx::SelectionList sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").storeVector(&sel).multiValue()));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    ASSERT_NO_THROW_GMX(assigner.startOption("sf"));
    // TODO: Should this be changed to an InvalidInputError?
    EXPECT_THROW_GMX(assigner.appendValue("nonexistentfile"), gmx::FileIOError);
    EXPECT_THROW_GMX(assigner.appendValue(TestFileManager::getInputFilePath("selfile.dat").string()),
                     gmx::InvalidInputError);
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
}


TEST_F(SelectionFileOptionTest, GivesErrorWithMultipleFiles)
{
    gmx::SelectionList sel;
    using gmx::SelectionOption;
    ASSERT_NO_THROW_GMX(options_.addOption(SelectionOption("sel").storeVector(&sel).multiValue()));

    gmx::OptionsAssigner assigner(&options_);
    EXPECT_NO_THROW_GMX(assigner.start());
    ASSERT_NO_THROW_GMX(assigner.startOption("sel"));
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    ASSERT_NO_THROW_GMX(assigner.startOption("sf"));
    EXPECT_NO_THROW_GMX(assigner.appendValue(TestFileManager::getInputFilePath("selfile.dat").string()));
    EXPECT_THROW_GMX(assigner.appendValue("nonexistentfile"), gmx::InvalidInputError);
    EXPECT_NO_THROW_GMX(assigner.finishOption());
    EXPECT_NO_THROW_GMX(assigner.finish());
    EXPECT_NO_THROW_GMX(options_.finish());
}

} // namespace
} // namespace test
} // namespace gmx
