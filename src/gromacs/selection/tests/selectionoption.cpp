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
 * Tests selection parsing and compilation.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/errorreporting/emptyerrorreporter.h"
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
    _sc.init();
    _sc.setReferencePosType("atom");
    _sc.setOutputPosType("atom");
    _options.globalProperties().setSelectionCollection(&_sc);
}


TEST_F(SelectionOptionTest, ParsesSimpleSelection)
{
    gmx::Selection *sel = NULL;
    using gmx::SelectionOption;
    _options.addOption(SelectionOption("sel").store(&sel));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&_options, &errors);
    ASSERT_EQ(0, assigner.startOption("sel"));
    EXPECT_EQ(0, assigner.appendValue("resname RA RB"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, _options.finish(&errors));
    ASSERT_TRUE(sel != NULL);
    ASSERT_FALSE(sel->isDynamic());
}


TEST_F(SelectionOptionTest, HandlesDynamicSelectionWhenStaticRequired)
{
    gmx::Selection *sel = NULL;
    using gmx::SelectionOption;
    _options.addOption(SelectionOption("sel").store(&sel).onlyStatic());

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&_options, &errors);
    ASSERT_EQ(0, assigner.startOption("sel"));
    EXPECT_NE(0, assigner.appendValue("resname RA RB and x < 5"));
    EXPECT_NE(0, assigner.finish());
    EXPECT_EQ(0, _options.finish(&errors));
}


TEST_F(SelectionOptionTest, HandlesTooManySelections)
{
    gmx::Selection *sel = NULL;
    using gmx::SelectionOption;
    _options.addOption(SelectionOption("sel").store(&sel));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&_options, &errors);
    ASSERT_EQ(0, assigner.startOption("sel"));
    EXPECT_EQ(0, assigner.appendValue("resname RA RB"));
    EXPECT_NE(0, assigner.appendValue("resname RB RC"));
    EXPECT_NE(0, assigner.finish());
    EXPECT_EQ(0, _options.finish(&errors));
    ASSERT_TRUE(sel != NULL);
}


TEST_F(SelectionOptionTest, HandlesTooFewSelections)
{
    gmx::Selection *sel[2] = {NULL, NULL};
    using gmx::SelectionOption;
    _options.addOption(SelectionOption("sel").store(sel).valueCount(2));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&_options, &errors);
    ASSERT_EQ(0, assigner.startOption("sel"));
    EXPECT_EQ(0, assigner.appendValue("resname RA RB"));
    EXPECT_NE(0, assigner.finish());
    EXPECT_EQ(0, _options.finish(&errors));
}


TEST_F(SelectionOptionTest, HandlesDelayedRequiredSelection)
{
    gmx::Selection *sel = NULL;
    using gmx::SelectionOption;
    _options.addOption(SelectionOption("sel").store(&sel).required());

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&_options, &errors);
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, _options.finish(&errors));
    EXPECT_EQ(0, _sc.parseRequestedFromString("resname RA RB"));
    ASSERT_TRUE(sel != NULL);
}


TEST_F(SelectionOptionTest, HandlesTooFewDelayedRequiredSelections)
{
    gmx::Selection *sel[2] = {NULL, NULL};
    using gmx::SelectionOption;
    _options.addOption(SelectionOption("sel").store(sel).required()
                           .valueCount(2));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&_options, &errors);
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, _options.finish(&errors));
    EXPECT_NE(0, _sc.parseRequestedFromString("resname RA RB"));
}


TEST_F(SelectionOptionTest, HandlesDelayedOptionalSelection)
{
    gmx::Selection *sel = NULL;
    using gmx::SelectionOption;
    _options.addOption(SelectionOption("sel").store(&sel));

    gmx::EmptyErrorReporter errors;
    gmx::OptionsAssigner assigner(&_options, &errors);
    ASSERT_EQ(0, assigner.startOption("sel"));
    EXPECT_EQ(0, assigner.finish());
    EXPECT_EQ(0, _options.finish(&errors));
    EXPECT_EQ(0, _sc.parseRequestedFromString("resname RA RB"));
    ASSERT_TRUE(sel != NULL);
}

} // namespace
