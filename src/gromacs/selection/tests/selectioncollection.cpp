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

#include "smalloc.h"
#include "tpxio.h"

#include "gromacs/errorreporting/emptyerrorreporter.h"
#include "gromacs/selection/poscalc.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selection.h"

namespace
{

class SelectionCollectionTest : public ::testing::Test
{
    public:
        SelectionCollectionTest();
        ~SelectionCollectionTest();

        void setAtomCount(int natoms)
        {
            _sc.setTopology(NULL, natoms);
        }
        void loadTopology(const char *filename);

        gmx::SelectionCollection _sc;
        gmx::EmptyErrorReporter  _errors;
        t_topology              *_top;
};

SelectionCollectionTest::SelectionCollectionTest()
    : _sc(NULL), _top(NULL)
{
    _sc.init();
    _sc.setReferencePosType("atom");
    _sc.setOutputPosType("atom");
}


SelectionCollectionTest::~SelectionCollectionTest()
{
    if (_top != NULL)
    {
        done_top(_top);
        sfree(_top);
    }
}


void
SelectionCollectionTest::loadTopology(const char *filename)
{
    char    title[STRLEN];
    int     ePBC;
    rvec   *xtop;
    matrix  box;

    snew(_top, 1);
    read_tps_conf(filename, title, _top, &ePBC, &xtop, NULL, box, FALSE);
    sfree(xtop);

    ASSERT_EQ(0, _sc.setTopology(_top, -1));
}


TEST_F(SelectionCollectionTest, HandlesNoSelections)
{
    EXPECT_FALSE(_sc.requiresTopology());
    EXPECT_EQ(0, _sc.compile());
}


TEST_F(SelectionCollectionTest, ParsesSimpleSelections)
{
    std::vector<gmx::Selection *> sel;
    EXPECT_EQ(0, _sc.parseFromString("atomnr 1 to 5", &_errors, &sel));
    EXPECT_EQ(1U, sel.size());
    EXPECT_EQ(0, _sc.parseFromString("atomnr 2 to 4", &_errors, &sel));
    EXPECT_EQ(2U, sel.size());
    setAtomCount(10);
    EXPECT_EQ(0, _sc.compile());
    EXPECT_EQ(5, sel[0]->posCount());
    EXPECT_EQ(3, sel[1]->posCount());
}


TEST_F(SelectionCollectionTest, ParsesArithmeticExpressions)
{
    std::vector<gmx::Selection *> sel;
    EXPECT_EQ(0, _sc.parseFromString("x+1 > 3", &_errors, &sel));
    EXPECT_EQ(1U, sel.size());
    EXPECT_EQ(0, _sc.parseFromString("(y-1)^2 <= 1", &_errors, &sel));
    EXPECT_EQ(2U, sel.size());
    loadTopology(SOURCE_DIR "/src/gromacs/selection/tests/simple.gro");
    EXPECT_EQ(0, _sc.compile());
    EXPECT_EQ(15, sel[0]->posCount());
    EXPECT_EQ(15, sel[1]->posCount());
}

} // namespace
