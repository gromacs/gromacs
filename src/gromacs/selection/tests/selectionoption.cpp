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
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/options/globalproperties.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoption.h"

namespace
{

class SelectionOptionTest : public ::testing::Test
{
    public:
        SelectionOptionTest();

        gmx::SelectionCollection _sc;
        gmx::Options             _options;
};

SelectionOptionTest::SelectionOptionTest()
    : _sc(NULL), _options(NULL, NULL)
{
    _sc.setReferencePosType("atom");
    _sc.setOutputPosType("atom");
    _options.globalProperties().setSelectionCollection(&_sc);
}


TEST_F(SelectionOptionTest, ParsesSimpleSelection)
{
    gmx::Selection *sel = NULL;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(_options.addOption(SelectionOption("sel").store(&sel)));

    gmx::OptionsAssigner assigner(&_options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("resname RA RB"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(_options.finish());

    ASSERT_TRUE(sel != NULL);
    ASSERT_FALSE(sel->isDynamic());
}


TEST_F(SelectionOptionTest, HandlesDynamicSelectionWhenStaticRequired)
{
    gmx::Selection *sel = NULL;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(_options.addOption(
                        SelectionOption("sel").store(&sel).onlyStatic()));

    gmx::OptionsAssigner assigner(&_options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_THROW(assigner.appendValue("resname RA RB and x < 5"), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(_options.finish());
}


TEST_F(SelectionOptionTest, HandlesTooManySelections)
{
    gmx::Selection *sel = NULL;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(_options.addOption(SelectionOption("sel").store(&sel)));

    gmx::OptionsAssigner assigner(&_options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("resname RA RB"));
    EXPECT_THROW(assigner.appendValue("resname RB RC"), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(_options.finish());
    ASSERT_TRUE(sel != NULL);
}


TEST_F(SelectionOptionTest, HandlesTooFewSelections)
{
    gmx::Selection *sel[2] = {NULL, NULL};
    using gmx::SelectionOption;
    ASSERT_NO_THROW(_options.addOption(
                        SelectionOption("sel").store(sel).valueCount(2)));

    gmx::OptionsAssigner assigner(&_options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("resname RA RB"));
    EXPECT_THROW(assigner.finishOption(), gmx::InvalidInputError);
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(_options.finish());
}


TEST_F(SelectionOptionTest, HandlesAdjuster)
{
    std::vector<gmx::Selection *> sel;
    gmx::SelectionOptionAdjuster *adjuster;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(_options.addOption(
                        SelectionOption("sel").storeVector(&sel).multiValue()
                            .getAdjuster(&adjuster)));

    gmx::OptionsAssigner assigner(&_options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("resname RA RB"));
    EXPECT_NO_THROW(assigner.appendValue("resname RB RC"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(_options.finish());
    EXPECT_NO_THROW(adjuster->setValueCount(2));
}


TEST_F(SelectionOptionTest, HandlesDynamicWhenStaticRequiredWithAdjuster)
{
    gmx::Selection *sel;
    gmx::SelectionOptionAdjuster *adjuster;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(_options.addOption(
                        SelectionOption("sel").store(&sel)
                            .getAdjuster(&adjuster)));

    gmx::OptionsAssigner assigner(&_options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("x < 5"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(_options.finish());
    EXPECT_THROW(adjuster->setOnlyStatic(true), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, HandlesTooManySelectionsWithAdjuster)
{
    std::vector<gmx::Selection *> sel;
    gmx::SelectionOptionAdjuster *adjuster;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(_options.addOption(
                        SelectionOption("sel").storeVector(&sel).multiValue()
                            .getAdjuster(&adjuster)));

    gmx::OptionsAssigner assigner(&_options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("resname RA RB"));
    EXPECT_NO_THROW(assigner.appendValue("resname RB RC"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(_options.finish());
    EXPECT_THROW(adjuster->setValueCount(1), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, HandlesTooFewSelectionsWithAdjuster)
{
    std::vector<gmx::Selection *> sel;
    gmx::SelectionOptionAdjuster *adjuster;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(_options.addOption(
                        SelectionOption("sel").storeVector(&sel).multiValue()
                            .getAdjuster(&adjuster)));

    gmx::OptionsAssigner assigner(&_options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.appendValue("resname RA RB"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(_options.finish());
    EXPECT_THROW(adjuster->setValueCount(2), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, HandlesDelayedRequiredSelection)
{
    gmx::Selection *sel = NULL;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(_options.addOption(
                        SelectionOption("sel").store(&sel).required()));

    gmx::OptionsAssigner assigner(&_options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(_options.finish());
    EXPECT_NO_THROW(_sc.parseRequestedFromString("resname RA RB"));
    ASSERT_TRUE(sel != NULL);
}


TEST_F(SelectionOptionTest, HandlesTooFewDelayedRequiredSelections)
{
    gmx::Selection *sel[2] = {NULL, NULL};
    using gmx::SelectionOption;
    ASSERT_NO_THROW(_options.addOption(
                        SelectionOption("sel").store(sel).required()
                            .valueCount(2)));

    gmx::OptionsAssigner assigner(&_options);
    EXPECT_NO_THROW(assigner.start());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(_options.finish());
    EXPECT_THROW(_sc.parseRequestedFromString("resname RA RB"), gmx::InvalidInputError);
}


TEST_F(SelectionOptionTest, HandlesDelayedOptionalSelection)
{
    gmx::Selection *sel = NULL;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(_options.addOption(SelectionOption("sel").store(&sel)));

    gmx::OptionsAssigner assigner(&_options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(_options.finish());
    EXPECT_NO_THROW(_sc.parseRequestedFromString("resname RA RB"));
    ASSERT_TRUE(sel != NULL);
}


TEST_F(SelectionOptionTest, HandlesDelayedSelectionWithAdjuster)
{
    std::vector<gmx::Selection *> sel;
    gmx::SelectionOptionAdjuster *adjuster;
    using gmx::SelectionOption;
    ASSERT_NO_THROW(_options.addOption(
                        SelectionOption("sel").storeVector(&sel).valueCount(3)
                            .getAdjuster(&adjuster)));

    gmx::OptionsAssigner assigner(&_options);
    EXPECT_NO_THROW(assigner.start());
    ASSERT_NO_THROW(assigner.startOption("sel"));
    EXPECT_NO_THROW(assigner.finishOption());
    EXPECT_NO_THROW(assigner.finish());
    EXPECT_NO_THROW(_options.finish());
    EXPECT_NO_THROW(adjuster->setValueCount(2));
    EXPECT_NO_THROW(_sc.parseRequestedFromString("resname RA RB; resname RB RC"));
}

} // namespace
