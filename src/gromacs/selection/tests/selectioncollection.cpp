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
#include "statutil.h"
#include "tpxio.h"
#include "vec.h"

#include "gromacs/errorreporting/emptyerrorreporter.h"
#include "gromacs/selection/poscalc.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selection.h"
#include "gromacs/utility/flags.h"

#include "testutils/datapath.h"
#include "testutils/refdata.h"

namespace
{

/********************************************************************
 * Test fixture for selection testing
 */

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
        t_trxframe              *_frame;
};


SelectionCollectionTest::SelectionCollectionTest()
    : _sc(NULL), _top(NULL), _frame(NULL)
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

    if (_frame != NULL)
    {
        sfree(_frame->x);
        sfree(_frame);
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
    read_tps_conf(gmx::test::getTestFilePath(filename).c_str(),
                  title, _top, &ePBC, &xtop, NULL, box, FALSE);

    snew(_frame, 1);
    _frame->flags  = TRX_NEED_X;
    _frame->natoms = _top->atoms.nr;
    _frame->bX     = TRUE;
    snew(_frame->x, _frame->natoms);
    memcpy(_frame->x, xtop, sizeof(*_frame->x) * _frame->natoms);
    _frame->bBox   = TRUE;
    copy_mat(box, _frame->box);

    ASSERT_EQ(0, _sc.setTopology(_top, -1));
}


/********************************************************************
 * Test fixture for selection testing with reference data
 */

class SelectionCollectionDataTest : public SelectionCollectionTest
{
    public:
        enum TestFlag
        {
            efTestPositionBlocks = 1<<0,
            efTestEvaluation     = 1<<1,
            efTestPositions      = 1<<2
        };
        typedef gmx::FlagsTemplate<TestFlag> TestFlags;

        SelectionCollectionDataTest() : _count(0), _framenr(0) {}

        void setFlags(TestFlags flags) { _flags = flags; }

        void runTest(int natoms, const char *const *selections);
        void runTest(const char *filename, const char *const *selections);

    private:
        void runParser(const char *const *selections);
        void runCompiler();
        void checkCompiled();
        void runEvaluate();
        void runEvaluateFinal();

        gmx::test::TestReferenceData  _data;
        std::vector<gmx::Selection *> _sel;
        size_t                        _count;
        int                           _framenr;
        TestFlags                     _flags;
};


void
SelectionCollectionDataTest::runParser(const char *const *selections)
{
    size_t varcount = 0;
    _count = 0;
    for (size_t i = 0; selections[i] != NULL; ++i)
    {
        ASSERT_EQ(0, _sc.parseFromString(selections[i], &_errors, &_sel));
        char buf[50];
        if (_sel.size() == _count)
        {
            snprintf(buf, 50, "Variable%dParse", static_cast<int>(varcount + 1));
            _data.startCompound("VariableParse", buf);
            _data.checkString(selections[i], "Input");
            _data.finishCompound();
            ++varcount;
        }
        else
        {
            snprintf(buf, 50, "Selection%dParse", static_cast<int>(_count + 1));
            _data.startCompound("SelectionParse", buf);
            _data.checkString(selections[i], "Input");
            _data.checkString(_sel[_count]->name(), "Name");
            _data.checkString(_sel[_count]->selectionText(), "Text");
            _data.checkBoolean(_sel[_count]->isDynamic(), "Dynamic");
            _data.finishCompound();
            ++_count;
        }
    }
}


void
SelectionCollectionDataTest::runCompiler()
{
    ASSERT_EQ(0, _sc.compile());
    ASSERT_EQ(_count, _sel.size());
    checkCompiled();
}


void
SelectionCollectionDataTest::checkCompiled()
{
    for (size_t i = 0; i < _count; ++i)
    {
        char buf[50];
        snprintf(buf, 50, "Selection%dCompile", static_cast<int>(i + 1));
        _data.startCompound("SelectionCompile", buf);
        if (_sel[i]->indexGroup() != NULL)
        {
            _data.checkSequenceInteger(_sel[i]->indexGroup()->isize,
                                       _sel[i]->indexGroup()->index,
                                       "Atoms");
        }
        else
        {
            _data.checkSequenceInteger(0, NULL, "Atoms");
        }
        if (_flags.test(efTestPositionBlocks))
        {
            _data.checkSequenceInteger(_sel[i]->posCount() + 1,
                                       _sel[i]->positions()->m.mapb.index,
                                       "PositionBlocks");
        }
        _data.finishCompound();
    }
}


void
SelectionCollectionDataTest::runEvaluate()
{
    ++_framenr;
    ASSERT_EQ(0, _sc.evaluate(_frame, NULL));
    for (size_t i = 0; i < _count; ++i)
    {
        char buf[50];
        snprintf(buf, 50, "Selection%dFrame%d",
                 static_cast<int>(i + 1), _framenr);
        _data.startCompound("SelectionFrame", buf);
        if (_sel[i]->indexGroup() != NULL)
        {
            _data.checkSequenceInteger(_sel[i]->indexGroup()->isize,
                                       _sel[i]->indexGroup()->index,
                                       "Atoms");
        }
        else
        {
            _data.checkSequenceInteger(0, NULL, "Atoms");
        }
        if (_flags.test(efTestPositionBlocks))
        {
            _data.checkSequenceInteger(_sel[i]->posCount() + 1,
                                       _sel[i]->positions()->m.mapb.index,
                                       "PositionBlocks");
        }
        if (_flags.test(efTestPositions))
        {
            _data.checkSequenceVector(_sel[i]->posCount(),
                                      _sel[i]->positions()->x,
                                       "Positions");
        }
        _data.finishCompound();
    }
}


void
SelectionCollectionDataTest::runEvaluateFinal()
{
    ASSERT_EQ(0, _sc.evaluateFinal(_framenr));
    if (!_data.isWriteMode())
    {
        checkCompiled();
    }
}


void
SelectionCollectionDataTest::runTest(int natoms, const char * const *selections)
{
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    setAtomCount(natoms);
    ASSERT_NO_FATAL_FAILURE(runCompiler());
}


void
SelectionCollectionDataTest::runTest(const char *filename, const char * const *selections)
{
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    loadTopology(filename);
    ASSERT_NO_FATAL_FAILURE(runCompiler());
    if (_flags.test(efTestEvaluation))
    {
        ASSERT_NO_FATAL_FAILURE(runEvaluate());
        ASSERT_NO_FATAL_FAILURE(runEvaluateFinal());
    }
}


/********************************************************************
 * Tests for SelectionCollection functionality without reference data
 */

TEST_F(SelectionCollectionTest, HandlesNoSelections)
{
    EXPECT_FALSE(_sc.requiresTopology());
    EXPECT_EQ(0, _sc.compile());
}


/********************************************************************
 * Tests for selection keywords
 */

TEST_F(SelectionCollectionDataTest, HandlesAllNone)
{
    static const char * const selections[] = {
        "all",
        "none",
        NULL
    };
    runTest(10, selections);
}

TEST_F(SelectionCollectionDataTest, HandlesAtomnr)
{
    static const char * const selections[] = {
        "atomnr 1 to 3 6 to 8",
        "atomnr 4 2 5 to 7",
        "atomnr <= 5",
        NULL
    };
    runTest(10, selections);
}

TEST_F(SelectionCollectionDataTest, HandlesResnr)
{
    static const char * const selections[] = {
        "resnr 1 2 5",
        "resnr 4 to 3",
        NULL
    };
    runTest("simple.gro", selections);
}

// TODO: Add test for "resindex"
// TODO: Add test for "molindex"

TEST_F(SelectionCollectionDataTest, HandlesAtomname)
{
    static const char * const selections[] = {
        "name CB",
        "name S1 S2",
        NULL
    };
    runTest("simple.gro", selections);
}

// TODO: Add test for atomtype
// TODO: Add test for insertcode
// TODO: Add test for chain
// TODO: Add test for mass
// TODO: Add test for charge
// TODO: Add test for altloc
// TODO: Add test for occupancy
// TODO: Add test for beta

TEST_F(SelectionCollectionDataTest, HandlesResname)
{
    static const char * const selections[] = {
        "resname RA",
        "resname RB RC",
        NULL
    };
    runTest("simple.gro", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesCoordinateKeywords)
{
    static const char * const selections[] = {
        "x < 3",
        "y >= 3",
        NULL
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositions);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesSameResidue)
{
    static const char * const selections[] = {
        "same residue as atomnr 1 4 12",
        NULL
    };
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesSameResidueName)
{
    static const char * const selections[] = {
        "same resname as atomnr 1 14",
        NULL
    };
    runTest("simple.gro", selections);
}


/********************************************************************
 * Tests for selection syntactic constructs
 */

TEST_F(SelectionCollectionDataTest, HandlesConstantPositions)
{
    static const char * const selections[] = {
        "[1, -2, 3.5]",
        NULL
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositions);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesWithinConstantPositions)
{
    static const char * const selections[] = {
        "within 1 of [2, 1, 0]",
        NULL
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositions);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesRegexMatching)
{
    static const char * const selections[] = {
        "resname \"R[BD]\"",
        NULL
    };
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesBasicBoolean)
{
    static const char * const selections[] = {
        "atomnr 1 to 5 and atomnr 2 to 7",
        "atomnr 1 to 5 or not atomnr 3 to 8",
        "atomnr 1 to 5 and atomnr 2 to 6 and not not atomnr 3 to 7",
        NULL
    };
    runTest(10, selections);
}


TEST_F(SelectionCollectionDataTest, HandlesArithmeticExpressions)
{
    static const char * const selections[] = {
        "x+1 > 3",
        "(y-1)^2 <= 1",
        "x+--1 > 3",
        "-x+-1 < -3",
        NULL
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositions);
    runTest("simple.gro", selections);
}


/********************************************************************
 * Tests for complex boolean syntax
 */

TEST_F(SelectionCollectionDataTest, HandlesBooleanStaticAnalysis)
{
    static const char * const selections[] = {
        "atomnr 1 to 5 and atomnr 2 to 7 and x < 2",
        "atomnr 1 to 5 and (atomnr 4 to 7 or x < 2)",
        "atomnr 1 to 5 and y < 3 and (atomnr 4 to 7 or x < 2)",
        "atomnr 1 to 5 and not (atomnr 4 to 7 or x < 2)",
        NULL
    };
    runTest(10, selections);
}


TEST_F(SelectionCollectionDataTest, HandlesBooleanStaticAnalysisWithVariables)
{
    static const char * const selections[] = {
        "foo = atomnr 4 to 7 or x < 2",
        "atomnr 1 to 4 and foo",
        "atomnr 2 to 6 and y < 3 and foo",
        "atomnr 6 to 10 and not foo",
        NULL
    };
    runTest(10, selections);
}


TEST_F(SelectionCollectionDataTest, HandlesBooleanStaticAnalysisWithMoreVariables)
{
    static const char * const selections[] = {
        "foo = atomnr 4 to 7",
        "bar = foo and x < 2",
        "bar2 = foo and y < 2",
        "atomnr 1 to 4 and bar",
        "atomnr 2 to 6 and y < 3 and bar2",
        "atomnr 6 to 10 and not foo",
        NULL
    };
    runTest(10, selections);
}

} // namespace
