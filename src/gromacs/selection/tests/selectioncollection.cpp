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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gtest/gtest.h>

#include "smalloc.h"
#include "statutil.h"
#include "tpxio.h"
#include "vec.h"

#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/poscalc.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selection.h"
#include "gromacs/utility/flags.h"
#include "gromacs/utility/format.h"

#include "testutils/datapath.h"
#include "testutils/refdata.h"
#include "testutils/testoptions.h"

namespace
{

/********************************************************************
 * Test fixture for selection testing
 */

class SelectionCollectionTest : public ::testing::Test
{
    public:
        static void SetUpTestCase();

        static int               s_debugLevel;

        SelectionCollectionTest();
        ~SelectionCollectionTest();

        void setAtomCount(int natoms)
        {
            ASSERT_NO_THROW(_sc.setTopology(NULL, natoms));
        }
        void loadTopology(const char *filename);

        gmx::SelectionCollection _sc;
        gmx::SelectionList       _sel;
        t_topology              *_top;
        t_trxframe              *_frame;
};

int SelectionCollectionTest::s_debugLevel = 0;

void SelectionCollectionTest::SetUpTestCase()
{
    gmx::Options options(NULL, NULL);
    options.addOption(gmx::IntegerOption("seldebug").store(&s_debugLevel));
    gmx::test::parseTestOptions(&options);
}


SelectionCollectionTest::SelectionCollectionTest()
    : _top(NULL), _frame(NULL)
{
    _sc.setDebugLevel(s_debugLevel);
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

    ASSERT_NO_THROW(_sc.setTopology(_top, -1));
}


/********************************************************************
 * Test fixture for selection testing with reference data
 */

class SelectionCollectionDataTest : public SelectionCollectionTest
{
    public:
        enum TestFlag
        {
            efTestEvaluation            = 1<<0,
            efTestPositionAtoms         = 1<<1,
            efTestPositionCoordinates   = 1<<2
        };
        typedef gmx::FlagsTemplate<TestFlag> TestFlags;

        SelectionCollectionDataTest()
            : _checker(_data.rootChecker()), _count(0), _framenr(0)
        {
        }

        void setFlags(TestFlags flags) { _flags = flags; }

        void runTest(int natoms, const char *const *selections);
        void runTest(const char *filename, const char *const *selections);

    private:
        static void checkSelection(gmx::test::TestReferenceChecker *checker,
                                   const gmx::Selection &sel, TestFlags flags);

        void runParser(const char *const *selections);
        void runCompiler();
        void checkCompiled();
        void runEvaluate();
        void runEvaluateFinal();

        gmx::test::TestReferenceData  _data;
        gmx::test::TestReferenceChecker _checker;
        size_t                        _count;
        int                           _framenr;
        TestFlags                     _flags;
};


void
SelectionCollectionDataTest::checkSelection(
        gmx::test::TestReferenceChecker *checker,
        const gmx::Selection &sel, TestFlags flags)
{
    using gmx::test::TestReferenceChecker;

    {
        gmx::ConstArrayRef<int> atoms = sel.atomIndices();
        checker->checkSequence(atoms.begin(), atoms.end(), "Atoms");
    }
    if (flags.test(efTestPositionAtoms)
        || flags.test(efTestPositionCoordinates))
    {
        TestReferenceChecker compound(
                checker->checkSequenceCompound("Positions", sel.posCount()));
        for (int i = 0; i < sel.posCount(); ++i)
        {
            TestReferenceChecker poscompound(compound.checkCompound("Position", NULL));
            const gmx::SelectionPosition &p = sel.position(i);
            if (flags.test(efTestPositionAtoms))
            {
                gmx::ConstArrayRef<int> atoms = p.atomIndices();
                poscompound.checkSequence(atoms.begin(), atoms.end(), "Atoms");
            }
            if (flags.test(efTestPositionCoordinates))
            {
                poscompound.checkVector(p.x(), "Coordinates");
            }
        }
    }
}


void
SelectionCollectionDataTest::runParser(const char *const *selections)
{
    using gmx::test::TestReferenceChecker;

    TestReferenceChecker compound(_checker.checkCompound("ParsedSelections", "Parsed"));
    size_t varcount = 0;
    _count = 0;
    for (size_t i = 0; selections[i] != NULL; ++i)
    {
        SCOPED_TRACE(std::string("Parsing selection \"")
                     + selections[i] + "\"");
        ASSERT_NO_THROW(_sc.parseFromString(selections[i], &_sel));
        if (_sel.size() == _count)
        {
            std::string id = gmx::formatString("Variable%d", static_cast<int>(varcount + 1));
            TestReferenceChecker varcompound(
                    compound.checkCompound("ParsedVariable", id.c_str()));
            varcompound.checkString(selections[i], "Input");
            ++varcount;
        }
        else
        {
            std::string id = gmx::formatString("Selection%d", static_cast<int>(_count + 1));
            TestReferenceChecker selcompound(
                    compound.checkCompound("ParsedSelection", id.c_str()));
            selcompound.checkString(selections[i], "Input");
            selcompound.checkString(_sel[_count].name(), "Name");
            selcompound.checkString(_sel[_count].selectionText(), "Text");
            selcompound.checkBoolean(_sel[_count].isDynamic(), "Dynamic");
            ++_count;
        }
    }
}


void
SelectionCollectionDataTest::runCompiler()
{
    ASSERT_NO_THROW(_sc.compile());
    ASSERT_EQ(_count, _sel.size());
    checkCompiled();
}


void
SelectionCollectionDataTest::checkCompiled()
{
    using gmx::test::TestReferenceChecker;
    const TestFlags mask = ~TestFlags(efTestPositionCoordinates);

    TestReferenceChecker compound(_checker.checkCompound("CompiledSelections", "Compiled"));
    for (size_t i = 0; i < _count; ++i)
    {
        SCOPED_TRACE(std::string("Checking selection \"") +
                     _sel[i].selectionText() + "\"");
        std::string id = gmx::formatString("Selection%d", static_cast<int>(i + 1));
        TestReferenceChecker selcompound(
                compound.checkCompound("Selection", id.c_str()));
        checkSelection(&selcompound, _sel[i], _flags & mask);
    }
}


void
SelectionCollectionDataTest::runEvaluate()
{
    using gmx::test::TestReferenceChecker;

    ++_framenr;
    ASSERT_NO_THROW(_sc.evaluate(_frame, NULL));
    std::string frame = gmx::formatString("Frame%d", _framenr);
    TestReferenceChecker compound(
            _checker.checkCompound("EvaluatedSelections", frame.c_str()));
    for (size_t i = 0; i < _count; ++i)
    {
        SCOPED_TRACE(std::string("Checking selection \"") +
                     _sel[i].selectionText() + "\"");
        std::string id = gmx::formatString("Selection%d", static_cast<int>(i + 1));
        TestReferenceChecker selcompound(
                compound.checkCompound("Selection", id.c_str()));
        checkSelection(&selcompound, _sel[i], _flags);
    }
}


void
SelectionCollectionDataTest::runEvaluateFinal()
{
    ASSERT_NO_THROW(_sc.evaluateFinal(_framenr));
    if (!_checker.isWriteMode())
    {
        checkCompiled();
    }
}


void
SelectionCollectionDataTest::runTest(int natoms, const char * const *selections)
{
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    ASSERT_NO_FATAL_FAILURE(setAtomCount(natoms));
    ASSERT_NO_FATAL_FAILURE(runCompiler());
}


void
SelectionCollectionDataTest::runTest(const char *filename, const char * const *selections)
{
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    ASSERT_NO_FATAL_FAILURE(loadTopology(filename));
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
    EXPECT_NO_THROW(_sc.compile());
}

TEST_F(SelectionCollectionTest, HandlesMissingMethodParamValue)
{
    EXPECT_THROW(_sc.parseFromString("mindist from atomnr 1 cutoff", &_sel),
                 gmx::InvalidInputError);
}

TEST_F(SelectionCollectionTest, HandlesMissingMethodParamValue2)
{
    EXPECT_THROW(_sc.parseFromString("within 1 of", &_sel),
                 gmx::InvalidInputError);
}

TEST_F(SelectionCollectionTest, HandlesMissingMethodParamValue3)
{
    EXPECT_THROW(_sc.parseFromString("within of atomnr 1", &_sel),
                 gmx::InvalidInputError);
}

// TODO: Tests for more parser errors

TEST_F(SelectionCollectionTest, RecoversFromUnknownGroupReference)
{
    ASSERT_NO_THROW(_sc.parseFromString("group \"foo\"", &_sel));
    ASSERT_NO_FATAL_FAILURE(setAtomCount(10));
    EXPECT_THROW(_sc.setIndexGroups(NULL), gmx::InvalidInputError);
    EXPECT_THROW(_sc.compile(), gmx::APIError);
}

TEST_F(SelectionCollectionTest, RecoversFromMissingMoleculeInfo)
{
    ASSERT_NO_THROW(_sc.parseFromString("molindex 1 to 5", &_sel));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW(_sc.compile(), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromMissingAtomTypes)
{
    ASSERT_NO_THROW(_sc.parseFromString("type CA", &_sel));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW(_sc.compile(), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromMissingPDBInfo)
{
    ASSERT_NO_THROW(_sc.parseFromString("altloc A", &_sel));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW(_sc.compile(), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromInvalidPermutation)
{
    ASSERT_NO_THROW(_sc.parseFromString("all permute 1 1", &_sel));
    ASSERT_NO_FATAL_FAILURE(setAtomCount(10));
    EXPECT_THROW(_sc.compile(), gmx::InvalidInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromInvalidPermutation2)
{
    ASSERT_NO_THROW(_sc.parseFromString("all permute 3 2 1", &_sel));
    ASSERT_NO_FATAL_FAILURE(setAtomCount(10));
    EXPECT_THROW(_sc.compile(), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromInvalidPermutation3)
{
    ASSERT_NO_THROW(_sc.parseFromString("x < 1.5 permute 3 2 1", &_sel));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    ASSERT_NO_THROW(_sc.compile());
    EXPECT_THROW(_sc.evaluate(_frame, NULL), gmx::InconsistentInputError);
}

// TODO: Tests for evaluation errors


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
        "x {-1 to 2}",
        NULL
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
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


TEST_F(SelectionCollectionDataTest, HandlesPositionKeywords)
{
    static const char * const selections[] = {
        "cog of resnr 1 3",
        "res_cog of name CB and resnr 1 3",
        "whole_res_cog of name CB and resnr 1 3",
        "part_res_cog of x < 3",
        "dyn_res_cog of x < 3",
        NULL
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates
             | efTestPositionAtoms);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesDistanceKeyword)
{
    static const char * const selections[] = {
        "distance from cog of resnr 1 < 2",
        NULL
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesMinDistanceKeyword)
{
    static const char * const selections[] = {
        "mindistance from resnr 1 < 2",
        NULL
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesWithinKeyword)
{
    static const char * const selections[] = {
        "within 1 of resnr 2",
        NULL
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


// TODO: Add test for "insolidangle"


// TODO: Check the handling of mapped and reference IDs in the modifier tests
// below.

TEST_F(SelectionCollectionDataTest, HandlesPermuteModifier)
{
    static const char * const selections[] = {
        "all permute 3 1 2",
        "res_cog of resnr 1 to 4 permute 2 1",
        "name CB S1 and res_cog x < 3 permute 2 1",
        NULL
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates
             | efTestPositionAtoms);
    runTest("simple.gro", selections);
}


// TODO: Add tests for plus/merge on dynamic selections
// (can't remember whether it's actually implemented or not).

TEST_F(SelectionCollectionDataTest, HandlesPlusModifier)
{
    static const char * const selections[] = {
        "name S2 plus name S1",
        "res_cog of resnr 2 plus res_cog of resnr 1 plus res_cog of resnr 3",
        NULL
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates
             | efTestPositionAtoms);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesMergeModifier)
{
    static const char * const selections[] = {
        "name S2 merge name S1",
        "name S2 merge name S1 merge name CB",
        NULL
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
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
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesWithinConstantPositions)
{
    static const char * const selections[] = {
        "within 1 of [2, 1, 0]",
        NULL
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}

#ifdef HAVE_REGEX_H
TEST_F(SelectionCollectionDataTest, HandlesRegexMatching)
#else
TEST_F(SelectionCollectionDataTest, DISABLED_HandlesRegexMatching)
#endif
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
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
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
