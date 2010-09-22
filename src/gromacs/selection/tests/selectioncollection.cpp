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

#include "gromacs/selection/poscalc.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selection.h"

namespace
{

class SelectionCollectionTest : public ::testing::Test
{
    public:
        void SetUp();
        void TearDown();

        gmx_ana_poscalc_coll_t *_pcc;
        gmx::SelectionCollection *_sc;
};

void SelectionCollectionTest::SetUp()
{
    _sc = NULL;
    ASSERT_EQ(0, gmx_ana_poscalc_coll_create(&_pcc));
    ASSERT_EQ(0, gmx::SelectionCollection::create(&_sc, _pcc));
}

void SelectionCollectionTest::TearDown()
{
    delete _sc;
    gmx_ana_poscalc_coll_free(_pcc);
}

TEST_F(SelectionCollectionTest, ParsesSimpleSelections)
{
    std::vector<gmx::Selection *> sel;
    _sc->setReferencePosType("atom");
    _sc->setOutputPosType("atom");
    EXPECT_EQ(0, _sc->parseFromString("resname RA RB", &sel));
    EXPECT_EQ(1U, sel.size());
    EXPECT_EQ(0, _sc->parseFromString("name CA SB", &sel));
    EXPECT_EQ(2U, sel.size());
}

} // namespace
