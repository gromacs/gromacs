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
#include "config.h"
#endif

#include <gtest/gtest.h>

#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/statutil.h"
#include "gromacs/legacyheaders/tpxio.h"
#include "gromacs/legacyheaders/vec.h"

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selection.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/flags.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"
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
            ASSERT_NO_THROW(sc_.setTopology(NULL, natoms));
        }
        void loadTopology(const char *filename);

        gmx::SelectionCollection sc_;
        gmx::SelectionList       sel_;
        t_topology              *top_;
        t_trxframe              *frame_;
};

int SelectionCollectionTest::s_debugLevel = 0;

void SelectionCollectionTest::SetUpTestCase()
{
    gmx::Options options(NULL, NULL);
    options.addOption(gmx::IntegerOption("seldebug").store(&s_debugLevel));
    gmx::test::parseTestOptions(&options);
}


SelectionCollectionTest::SelectionCollectionTest()
    : top_(NULL), frame_(NULL)
{
    sc_.setDebugLevel(s_debugLevel);
    sc_.setReferencePosType("atom");
    sc_.setOutputPosType("atom");
}


SelectionCollectionTest::~SelectionCollectionTest()
{
    if (top_ != NULL)
    {
        done_top(top_);
        sfree(top_);
    }

    if (frame_ != NULL)
    {
        sfree(frame_->x);
        sfree(frame_);
    }
}


void
SelectionCollectionTest::loadTopology(const char *filename)
{
    char    title[STRLEN];
    int     ePBC;
    rvec   *xtop;
    matrix  box;

    snew(top_, 1);
    read_tps_conf(gmx::test::TestFileManager::getInputFilePath(filename).c_str(),
                  title, top_, &ePBC, &xtop, NULL, box, FALSE);

    snew(frame_, 1);
    frame_->flags  = TRX_NEED_X;
    frame_->natoms = top_->atoms.nr;
    frame_->bX     = TRUE;
    snew(frame_->x, frame_->natoms);
    memcpy(frame_->x, xtop, sizeof(*frame_->x) * frame_->natoms);
    frame_->bBox   = TRUE;
    copy_mat(box, frame_->box);

    ASSERT_NO_THROW(sc_.setTopology(top_, -1));
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
            : checker_(data_.rootChecker()), count_(0), framenr_(0)
        {
        }

        void setFlags(TestFlags flags) { flags_ = flags; }

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

        gmx::test::TestReferenceData  data_;
        gmx::test::TestReferenceChecker checker_;
        size_t                        count_;
        int                           framenr_;
        TestFlags                     flags_;
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

    TestReferenceChecker compound(checker_.checkCompound("ParsedSelections", "Parsed"));
    size_t varcount = 0;
    count_ = 0;
    for (size_t i = 0; selections[i] != NULL; ++i)
    {
        SCOPED_TRACE(std::string("Parsing selection \"")
                     + selections[i] + "\"");
        gmx::SelectionList result;
        ASSERT_NO_THROW(result = sc_.parseFromString(selections[i]));
        sel_.insert(sel_.end(), result.begin(), result.end());
        if (sel_.size() == count_)
        {
            std::string id = gmx::formatString("Variable%d", static_cast<int>(varcount + 1));
            TestReferenceChecker varcompound(
                    compound.checkCompound("ParsedVariable", id.c_str()));
            varcompound.checkString(selections[i], "Input");
            ++varcount;
        }
        else
        {
            std::string id = gmx::formatString("Selection%d", static_cast<int>(count_ + 1));
            TestReferenceChecker selcompound(
                    compound.checkCompound("ParsedSelection", id.c_str()));
            selcompound.checkString(selections[i], "Input");
            selcompound.checkString(sel_[count_].name(), "Name");
            selcompound.checkString(sel_[count_].selectionText(), "Text");
            selcompound.checkBoolean(sel_[count_].isDynamic(), "Dynamic");
            ++count_;
        }
    }
}


void
SelectionCollectionDataTest::runCompiler()
{
    ASSERT_NO_THROW(sc_.compile());
    ASSERT_EQ(count_, sel_.size());
    checkCompiled();
}


void
SelectionCollectionDataTest::checkCompiled()
{
    using gmx::test::TestReferenceChecker;
    const TestFlags mask = ~TestFlags(efTestPositionCoordinates);

    TestReferenceChecker compound(checker_.checkCompound("CompiledSelections", "Compiled"));
    for (size_t i = 0; i < count_; ++i)
    {
        SCOPED_TRACE(std::string("Checking selection \"") +
                     sel_[i].selectionText() + "\"");
        std::string id = gmx::formatString("Selection%d", static_cast<int>(i + 1));
        TestReferenceChecker selcompound(
                compound.checkCompound("Selection", id.c_str()));
        checkSelection(&selcompound, sel_[i], flags_ & mask);
    }
}


void
SelectionCollectionDataTest::runEvaluate()
{
    using gmx::test::TestReferenceChecker;

    ++framenr_;
    ASSERT_NO_THROW(sc_.evaluate(frame_, NULL));
    std::string frame = gmx::formatString("Frame%d", framenr_);
    TestReferenceChecker compound(
            checker_.checkCompound("EvaluatedSelections", frame.c_str()));
    for (size_t i = 0; i < count_; ++i)
    {
        SCOPED_TRACE(std::string("Checking selection \"") +
                     sel_[i].selectionText() + "\"");
        std::string id = gmx::formatString("Selection%d", static_cast<int>(i + 1));
        TestReferenceChecker selcompound(
                compound.checkCompound("Selection", id.c_str()));
        checkSelection(&selcompound, sel_[i], flags_);
    }
}


void
SelectionCollectionDataTest::runEvaluateFinal()
{
    ASSERT_NO_THROW(sc_.evaluateFinal(framenr_));
    if (!checker_.isWriteMode())
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
    if (flags_.test(efTestEvaluation))
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
    EXPECT_FALSE(sc_.requiresTopology());
    EXPECT_NO_THROW(sc_.compile());
}

TEST_F(SelectionCollectionTest, ParsesSelectionsFromFile)
{
    ASSERT_NO_THROW(sel_ = sc_.parseFromFile(
                gmx::test::TestFileManager::getInputFilePath("selfile.dat")));
    // These should match the contents of selfile.dat
    ASSERT_EQ(2U, sel_.size());
    EXPECT_STREQ("resname RA RB", sel_[0].selectionText());
    EXPECT_STREQ("resname RB RC", sel_[1].selectionText());
}

#ifdef HAVE_REGEX_H
TEST_F(SelectionCollectionTest, HandlesInvalidRegularExpressions)
{
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW({
            sc_.parseFromString("resname ~ \"R[A\"");
            sc_.compile();
        }, gmx::InvalidInputError);
}
#else
TEST_F(SelectionCollectionTest, HandlesUnsupportedRegularExpressions)
{
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW({
            sc_.parseFromString("resname \"R[AD]\"");
            sc_.compile();
        }, gmx::InvalidInputError);
}
#endif

TEST_F(SelectionCollectionTest, HandlesMissingMethodParamValue)
{
    EXPECT_THROW(sc_.parseFromString("mindist from atomnr 1 cutoff"),
                 gmx::InvalidInputError);
}

TEST_F(SelectionCollectionTest, HandlesMissingMethodParamValue2)
{
    EXPECT_THROW(sc_.parseFromString("within 1 of"),
                 gmx::InvalidInputError);
}

TEST_F(SelectionCollectionTest, HandlesMissingMethodParamValue3)
{
    EXPECT_THROW(sc_.parseFromString("within of atomnr 1"),
                 gmx::InvalidInputError);
}

TEST_F(SelectionCollectionTest, HandlesHelpKeywordInInvalidContext)
{
    EXPECT_THROW(sc_.parseFromString("resname help"),
                 gmx::InvalidInputError);
}

// TODO: Tests for more parser errors

TEST_F(SelectionCollectionTest, RecoversFromUnknownGroupReference)
{
    ASSERT_NO_THROW(sc_.parseFromString("group \"foo\""));
    ASSERT_NO_FATAL_FAILURE(setAtomCount(10));
    EXPECT_THROW(sc_.setIndexGroups(NULL), gmx::InvalidInputError);
    EXPECT_THROW(sc_.compile(), gmx::APIError);
}

TEST_F(SelectionCollectionTest, RecoversFromMissingMoleculeInfo)
{
    ASSERT_NO_THROW(sc_.parseFromString("molindex 1 to 5"));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW(sc_.compile(), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromMissingAtomTypes)
{
    ASSERT_NO_THROW(sc_.parseFromString("type CA"));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW(sc_.compile(), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromMissingPDBInfo)
{
    ASSERT_NO_THROW(sc_.parseFromString("altloc A"));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW(sc_.compile(), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromInvalidPermutation)
{
    ASSERT_NO_THROW(sc_.parseFromString("all permute 1 1"));
    ASSERT_NO_FATAL_FAILURE(setAtomCount(10));
    EXPECT_THROW(sc_.compile(), gmx::InvalidInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromInvalidPermutation2)
{
    ASSERT_NO_THROW(sc_.parseFromString("all permute 3 2 1"));
    ASSERT_NO_FATAL_FAILURE(setAtomCount(10));
    EXPECT_THROW(sc_.compile(), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromInvalidPermutation3)
{
    ASSERT_NO_THROW(sc_.parseFromString("x < 1.5 permute 3 2 1"));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    ASSERT_NO_THROW(sc_.compile());
    EXPECT_THROW(sc_.evaluate(frame_, NULL), gmx::InconsistentInputError);
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
        "resid 4 to 3",
        NULL
    };
    runTest("simple.gro", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesResIndex)
{
    static const char * const selections[] = {
        "resindex 1 4",
        "residue 1 3",
        NULL
    };
    runTest("simple.pdb", selections);
}

// TODO: Add test for "molindex"

TEST_F(SelectionCollectionDataTest, HandlesAtomname)
{
    static const char * const selections[] = {
        "name CB",
        "atomname S1 S2",
        NULL
    };
    runTest("simple.gro", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesPdbAtomname)
{
    static const char * const selections[] = {
        "name HG21",
        "name 1HG2",
        "pdbname HG21 CB",
        "pdbatomname 1HG2",
        NULL
    };
    runTest("simple.pdb", selections);
}

// TODO: Add test for atomtype

TEST_F(SelectionCollectionDataTest, HandlesChain)
{
    static const char * const selections[] = {
        "chain A",
        "chain B",
        NULL
    };
    runTest("simple.pdb", selections);
}

// TODO: Add test for mass
// TODO: Add test for charge

TEST_F(SelectionCollectionDataTest, HandlesAltLoc)
{
    static const char * const selections[] = {
        "altloc \" \"",
        "altloc A",
        NULL
    };
    runTest("simple.pdb", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesInsertCode)
{
    static const char * const selections[] = {
        "insertcode \" \"",
        "insertcode A",
        NULL
    };
    runTest("simple.pdb", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesOccupancy)
{
    static const char * const selections[] = {
        "occupancy 1",
        "occupancy < .5",
        NULL
    };
    runTest("simple.pdb", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesBeta)
{
    static const char * const selections[] = {
        "beta 0",
        "beta >= 0.3",
        NULL
    };
    runTest("simple.pdb", selections);
}

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


TEST_F(SelectionCollectionDataTest, HandlesForcedStringMatchingMode)
{
    static const char * const selections[] = {
        "name = S1 \"C?\"",
        "name ? S1 \"C?\"",
        NULL
    };
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesWildcardMatching)
{
    static const char * const selections[] = {
        "name \"S?\"",
        "name ? \"S?\"",
        NULL
    };
    runTest("simple.gro", selections);
}


#ifdef HAVE_REGEX_H
TEST_F(SelectionCollectionDataTest, HandlesRegexMatching)
{
    static const char * const selections[] = {
        "resname \"R[BD]\"",
        "resname ~ \"R[BD]\"",
        NULL
    };
    runTest("simple.gro", selections);
}
#endif


TEST_F(SelectionCollectionDataTest, HandlesBasicBoolean)
{
    static const char * const selections[] = {
        "atomnr 1 to 5 and atomnr 2 to 7",
        "atomnr 1 to 5 or not atomnr 3 to 8",
        "not not atomnr 1 to 5 and atomnr 2 to 6 and not not atomnr 3 to 7",
        "atomnr 1 to 5 and (atomnr 2 to 7 and atomnr 3 to 6)",
        "x < 5 and atomnr 1 to 5 and y < 3 and atomnr 2 to 4",
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


TEST_F(SelectionCollectionDataTest, HandlesNumericVariables)
{
    static const char * const selections[] = {
        "value = x + y",
        "value <= 4",
        "index = resnr",
        "index < 3",
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
        "atomnr 1 to 5 or (atomnr 4 to 6 and (atomnr 5 to 7 or x < 2))",
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


/********************************************************************
 * Tests for complex subexpression cases
 *
 * These tests use some knowledge of the implementation to trigger different
 * paths in the code.
 */

TEST_F(SelectionCollectionDataTest, HandlesUnusedVariables)
{
    static const char * const selections[] = {
        "unused1 = atomnr 1 to 3",
        "foo = atomnr 4 to 7",
        "atomnr 1 to 6 and foo",
        "unused2 = atomnr 3 to 5",
        NULL
    };
    runTest(10, selections);
}


TEST_F(SelectionCollectionDataTest, HandlesVariablesWithStaticEvaluationGroups)
{
    static const char * const selections[] = {
        "foo = atomnr 4 to 7 and x < 2",
        "atomnr 1 to 5 and foo",
        "atomnr 3 to 7 and foo",
        NULL
    };
    runTest(10, selections);
}


TEST_F(SelectionCollectionDataTest, HandlesVariablesWithMixedEvaluationGroups)
{
    static const char * const selections[] = {
        "foo = atomnr 4 to 7 and x < 2",
        "atomnr 1 to 6 and foo",
        "within 1 of foo",
        "foo",
        NULL
    };
    runTest(10, selections);
}


} // namespace
