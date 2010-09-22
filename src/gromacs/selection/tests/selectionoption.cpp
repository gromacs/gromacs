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
#include "gromacs/selection/poscalc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoption.h"

namespace
{

class SelectionOptionTest : public ::testing::Test
{
    public:
        SelectionOptionTest();
        ~SelectionOptionTest();

        gmx_ana_poscalc_coll_t *_pcc;
        gmx::SelectionCollection *_sc;
        gmx::Options            _options;
};

SelectionOptionTest::SelectionOptionTest()
    : _pcc(NULL), _sc(NULL), _options(NULL, NULL)
{
    gmx_ana_poscalc_coll_create(&_pcc);
    gmx::SelectionCollection::create(&_sc, _pcc);
    _sc->setReferencePosType("atom");
    _sc->setOutputPosType("atom");
    _options.globalProperties().setSelectionCollection(_sc);
}

SelectionOptionTest::~SelectionOptionTest()
{
    delete _sc;
    gmx_ana_poscalc_coll_free(_pcc);
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
}

} // namespace
