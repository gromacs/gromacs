/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Tests selection parsing and compilation.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/selection/selectioncollection.h"

#include <gtest/gtest.h>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/selection.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/flags.h"
#include "gromacs/utility/gmxregex.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/interactivetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testoptions.h"

#include "toputils.h"

namespace
{

/********************************************************************
 * Test fixture for selection testing
 */

class SelectionCollectionTest : public ::testing::Test
{
    public:
        static int               s_debugLevel;

        SelectionCollectionTest();
        ~SelectionCollectionTest();

        void setAtomCount(int natoms)
        {
            ASSERT_NO_THROW_GMX(sc_.setTopology(nullptr, natoms));
        }
        void loadTopology(const char *filename);
        void setTopology();
        void loadIndexGroups(const char *filename);

        gmx::test::TopologyManager  topManager_;
        gmx::SelectionCollection    sc_;
        gmx::SelectionList          sel_;
        gmx_ana_indexgrps_t        *grps_;
};

int SelectionCollectionTest::s_debugLevel = 0;

// cond/endcond do not seem to work here with Doxygen 1.8.5 parser.
#ifndef DOXYGEN
GMX_TEST_OPTIONS(SelectionCollectionTestOptions, options)
{
    options->addOption(gmx::IntegerOption("seldebug")
                           .store(&SelectionCollectionTest::s_debugLevel)
                           .description("Set selection debug level"));
}
#endif

SelectionCollectionTest::SelectionCollectionTest()
    : grps_(nullptr)
{
    topManager_.requestFrame();
    sc_.setDebugLevel(s_debugLevel);
    sc_.setReferencePosType("atom");
    sc_.setOutputPosType("atom");
}

SelectionCollectionTest::~SelectionCollectionTest()
{
    if (grps_ != nullptr)
    {
        gmx_ana_indexgrps_free(grps_);
    }
}

void
SelectionCollectionTest::loadTopology(const char *filename)
{
    topManager_.loadTopology(filename);
    setTopology();
}

void
SelectionCollectionTest::setTopology()
{
    ASSERT_NO_THROW_GMX(sc_.setTopology(topManager_.topology(), -1));
}

void
SelectionCollectionTest::loadIndexGroups(const char *filename)
{
    GMX_RELEASE_ASSERT(grps_ == nullptr,
                       "External groups can only be loaded once");
    std::string fullpath =
        gmx::test::TestFileManager::getInputFilePath(filename);
    gmx_ana_indexgrps_init(&grps_, nullptr, fullpath.c_str());
    sc_.setIndexGroups(grps_);
}


/********************************************************************
 * Test fixture for interactive SelectionCollection tests
 */

class SelectionCollectionInteractiveTest : public SelectionCollectionTest
{
    public:
        SelectionCollectionInteractiveTest()
            : helper_(data_.rootChecker())
        {
        }

        void runTest(int count, bool bInteractive,
                     const gmx::ArrayRef<const char *const> &input);

        gmx::test::TestReferenceData      data_;
        gmx::test::InteractiveTestHelper  helper_;
};

void SelectionCollectionInteractiveTest::runTest(
        int count, bool bInteractive,
        const gmx::ArrayRef<const char *const> &inputLines)
{
    helper_.setInputLines(inputLines);
    // TODO: Check something about the returned selections as well.
    ASSERT_NO_THROW_GMX(sc_.parseInteractive(
                                count, &helper_.inputStream(),
                                bInteractive ? &helper_.outputStream() : nullptr,
                                "for test context"));
    helper_.checkSession();
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
            efTestPositionCoordinates   = 1<<2,
            efTestPositionMapping       = 1<<3,
            efTestPositionMasses        = 1<<4,
            efTestPositionCharges       = 1<<5,
            efTestSelectionNames        = 1<<6,
            efDontTestCompiledAtoms     = 1<<8
        };
        typedef gmx::FlagsTemplate<TestFlag> TestFlags;

        SelectionCollectionDataTest()
            : checker_(data_.rootChecker()), count_(0), framenr_(0)
        {
        }

        void setFlags(TestFlags flags) { flags_ = flags; }

        void runParser(const gmx::ArrayRef<const char *const> &selections);
        void runCompiler();
        void runEvaluate();
        void runEvaluateFinal();

        void runTest(int                                     natoms,
                     const gmx::ArrayRef<const char *const> &selections);
        void runTest(const char                             *filename,
                     const gmx::ArrayRef<const char *const> &selections);

    private:
        static void checkSelection(gmx::test::TestReferenceChecker *checker,
                                   const gmx::Selection &sel, TestFlags flags);

        void checkCompiled();

        gmx::test::TestReferenceData    data_;
        gmx::test::TestReferenceChecker checker_;
        size_t                          count_;
        int                             framenr_;
        TestFlags                       flags_;
};


void
SelectionCollectionDataTest::checkSelection(
        gmx::test::TestReferenceChecker *checker,
        const gmx::Selection &sel, TestFlags flags)
{
    using gmx::test::TestReferenceChecker;

    {
        gmx::ArrayRef<const int> atoms = sel.atomIndices();
        checker->checkSequence(atoms.begin(), atoms.end(), "Atoms");
    }
    if (flags.test(efTestPositionAtoms)
        || flags.test(efTestPositionCoordinates)
        || flags.test(efTestPositionMapping)
        || flags.test(efTestPositionMasses)
        || flags.test(efTestPositionCharges))
    {
        TestReferenceChecker compound(
                checker->checkSequenceCompound("Positions", sel.posCount()));
        for (int i = 0; i < sel.posCount(); ++i)
        {
            TestReferenceChecker          poscompound(compound.checkCompound("Position", nullptr));
            const gmx::SelectionPosition &p = sel.position(i);
            if (flags.test(efTestPositionAtoms))
            {
                gmx::ArrayRef<const int> atoms = p.atomIndices();
                poscompound.checkSequence(atoms.begin(), atoms.end(), "Atoms");
            }
            if (flags.test(efTestPositionCoordinates))
            {
                poscompound.checkVector(p.x(), "Coordinates");
            }
            if (flags.test(efTestPositionMapping))
            {
                poscompound.checkInteger(p.refId(), "RefId");
                poscompound.checkInteger(p.mappedId(), "MappedId");
            }
            if (flags.test(efTestPositionMasses))
            {
                poscompound.checkReal(p.mass(), "Mass");
            }
            if (flags.test(efTestPositionCharges))
            {
                poscompound.checkReal(p.charge(), "Charge");
            }
        }
    }
}


void
SelectionCollectionDataTest::runParser(
        const gmx::ArrayRef<const char *const> &selections)
{
    using gmx::test::TestReferenceChecker;

    TestReferenceChecker compound(checker_.checkCompound("ParsedSelections", "Parsed"));
    size_t               varcount = 0;
    count_ = 0;
    for (size_t i = 0; i < selections.size(); ++i)
    {
        SCOPED_TRACE(std::string("Parsing selection \"")
                     + selections[i] + "\"");
        gmx::SelectionList result;
        ASSERT_NO_THROW_GMX(result = sc_.parseFromString(selections[i]));
        sel_.insert(sel_.end(), result.begin(), result.end());
        if (sel_.size() == count_)
        {
            std::string          id = gmx::formatString("Variable%d", static_cast<int>(varcount + 1));
            TestReferenceChecker varcompound(
                    compound.checkCompound("ParsedVariable", id.c_str()));
            varcompound.checkString(selections[i], "Input");
            ++varcount;
        }
        else
        {
            std::string          id = gmx::formatString("Selection%d", static_cast<int>(count_ + 1));
            TestReferenceChecker selcompound(
                    compound.checkCompound("ParsedSelection", id.c_str()));
            selcompound.checkString(selections[i], "Input");
            if (flags_.test(efTestSelectionNames))
            {
                selcompound.checkString(sel_[count_].name(), "Name");
            }
            selcompound.checkString(sel_[count_].selectionText(), "Text");
            selcompound.checkBoolean(sel_[count_].isDynamic(), "Dynamic");
            ++count_;
        }
    }
}


void
SelectionCollectionDataTest::runCompiler()
{
    ASSERT_NO_THROW_GMX(sc_.compile());
    ASSERT_EQ(count_, sel_.size());
    checkCompiled();
}


void
SelectionCollectionDataTest::checkCompiled()
{
    using gmx::test::TestReferenceChecker;
    const TestFlags      mask = ~TestFlags(efTestPositionCoordinates);

    TestReferenceChecker compound(checker_.checkCompound("CompiledSelections", "Compiled"));
    for (size_t i = 0; i < count_; ++i)
    {
        SCOPED_TRACE(std::string("Checking selection \"") +
                     sel_[i].selectionText() + "\"");
        std::string          id = gmx::formatString("Selection%d", static_cast<int>(i + 1));
        TestReferenceChecker selcompound(
                compound.checkCompound("Selection", id.c_str()));
        if (flags_.test(efTestSelectionNames))
        {
            selcompound.checkString(sel_[i].name(), "Name");
        }
        if (!flags_.test(efDontTestCompiledAtoms))
        {
            checkSelection(&selcompound, sel_[i], flags_ & mask);
        }
    }
}


void
SelectionCollectionDataTest::runEvaluate()
{
    using gmx::test::TestReferenceChecker;

    ++framenr_;
    ASSERT_NO_THROW_GMX(sc_.evaluate(topManager_.frame(), nullptr));
    std::string          frame = gmx::formatString("Frame%d", framenr_);
    TestReferenceChecker compound(
            checker_.checkCompound("EvaluatedSelections", frame.c_str()));
    for (size_t i = 0; i < count_; ++i)
    {
        SCOPED_TRACE(std::string("Checking selection \"") +
                     sel_[i].selectionText() + "\"");
        std::string          id = gmx::formatString("Selection%d", static_cast<int>(i + 1));
        TestReferenceChecker selcompound(
                compound.checkCompound("Selection", id.c_str()));
        checkSelection(&selcompound, sel_[i], flags_);
    }
}


void
SelectionCollectionDataTest::runEvaluateFinal()
{
    ASSERT_NO_THROW_GMX(sc_.evaluateFinal(framenr_));
    checkCompiled();
}


void
SelectionCollectionDataTest::runTest(
        int natoms, const gmx::ArrayRef<const char *const> &selections)
{
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    ASSERT_NO_FATAL_FAILURE(setAtomCount(natoms));
    ASSERT_NO_FATAL_FAILURE(runCompiler());
}


void
SelectionCollectionDataTest::runTest(
        const char *filename, const gmx::ArrayRef<const char *const> &selections)
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
    EXPECT_FALSE(sc_.requiredTopologyProperties().hasAny());
    EXPECT_NO_THROW_GMX(sc_.compile());
    EXPECT_FALSE(sc_.requiredTopologyProperties().hasAny());
}

TEST_F(SelectionCollectionTest, HandlesNoSelectionsWithDefaultPositionType)
{
    EXPECT_NO_THROW_GMX(sc_.setOutputPosType("res_com"));
    EXPECT_TRUE(sc_.requiredTopologyProperties().needsTopology);
    EXPECT_TRUE(sc_.requiredTopologyProperties().needsMasses);
    EXPECT_NO_THROW_GMX(sc_.setOutputPosType("res_cog"));
    EXPECT_TRUE(sc_.requiredTopologyProperties().needsTopology);
    EXPECT_FALSE(sc_.requiredTopologyProperties().needsMasses);
    ASSERT_NO_THROW_GMX(sc_.parseFromString("atom of atomnr 1 to 10"));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    ASSERT_NO_THROW_GMX(sc_.compile());
    EXPECT_FALSE(sc_.requiredTopologyProperties().hasAny());
}

TEST_F(SelectionCollectionTest, HandlesVelocityAndForceRequests)
{
    ASSERT_NO_THROW_GMX(sel_ = sc_.parseFromString("atomnr 1 to 10; none"));
    EXPECT_FALSE(sc_.requiredTopologyProperties().hasAny());
    ASSERT_NO_FATAL_FAILURE(setAtomCount(10));
    ASSERT_EQ(2U, sel_.size());
    ASSERT_NO_THROW_GMX(sel_[0].setEvaluateVelocities(true));
    ASSERT_NO_THROW_GMX(sel_[1].setEvaluateVelocities(true));
    ASSERT_NO_THROW_GMX(sel_[0].setEvaluateForces(true));
    ASSERT_NO_THROW_GMX(sel_[1].setEvaluateForces(true));
    EXPECT_FALSE(sc_.requiredTopologyProperties().hasAny());
    ASSERT_NO_THROW_GMX(sc_.compile());
    EXPECT_FALSE(sc_.requiredTopologyProperties().hasAny());
    EXPECT_TRUE(sel_[0].hasVelocities());
    EXPECT_TRUE(sel_[1].hasVelocities());
    EXPECT_TRUE(sel_[0].hasForces());
    EXPECT_TRUE(sel_[1].hasForces());
}

TEST_F(SelectionCollectionTest, HandlesForceRequestForCenterOfGeometry)
{
    ASSERT_NO_THROW_GMX(sel_ = sc_.parseFromString("res_cog of atomnr 1 to 10"));
    EXPECT_TRUE(sc_.requiredTopologyProperties().needsTopology);
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    ASSERT_EQ(1U, sel_.size());
    ASSERT_NO_THROW_GMX(sel_[0].setEvaluateForces(true));
    // In principle, the code could know here that the masses are required, but
    // currently it only knows this after compilation.
    ASSERT_NO_THROW_GMX(sc_.compile());
    EXPECT_TRUE(sc_.requiredTopologyProperties().needsMasses);
    EXPECT_TRUE(sel_[0].hasForces());
}

TEST_F(SelectionCollectionTest, ParsesSelectionsFromFile)
{
    ASSERT_NO_THROW_GMX(sel_ = sc_.parseFromFile(
                                    gmx::test::TestFileManager::getInputFilePath("selfile.dat")));
    // These should match the contents of selfile.dat
    ASSERT_EQ(2U, sel_.size());
    EXPECT_STREQ("resname RA RB", sel_[0].selectionText());
    EXPECT_STREQ("resname RB RC", sel_[1].selectionText());
}

TEST_F(SelectionCollectionTest, HandlesAtypicalWhitespace)
{
    ASSERT_NO_THROW_GMX(sel_ = sc_.parseFromString("atomnr\n1\r\nto\t10;\vatomnr 3\f to 14\r"));
    ASSERT_EQ(2U, sel_.size());
    EXPECT_STREQ("atomnr 1 to 10", sel_[0].selectionText());
    // TODO: Get rid of the trailing whitespace.
    EXPECT_STREQ("atomnr 3 to 14 ", sel_[1].selectionText());
}

TEST_F(SelectionCollectionTest, HandlesInvalidRegularExpressions)
{
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW_GMX({
                         sc_.parseFromString("resname ~ \"R[A\"");
                         sc_.compile();
                     }, gmx::InvalidInputError);
}

TEST_F(SelectionCollectionTest, HandlesUnsupportedRegularExpressions)
{
    if (!gmx::Regex::isSupported())
    {
        ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
        EXPECT_THROW_GMX({
                             sc_.parseFromString("resname \"R[AD]\"");
                             sc_.compile();
                         }, gmx::InvalidInputError);
    }
}

TEST_F(SelectionCollectionTest, HandlesMissingMethodParamValue)
{
    EXPECT_THROW_GMX(sc_.parseFromString("mindist from atomnr 1 cutoff"),
                     gmx::InvalidInputError);
}

TEST_F(SelectionCollectionTest, HandlesMissingMethodParamValue2)
{
    EXPECT_THROW_GMX(sc_.parseFromString("within 1 of"),
                     gmx::InvalidInputError);
}

TEST_F(SelectionCollectionTest, HandlesMissingMethodParamValue3)
{
    EXPECT_THROW_GMX(sc_.parseFromString("within of atomnr 1"),
                     gmx::InvalidInputError);
}

// TODO: Tests for more parser errors

TEST_F(SelectionCollectionTest, HandlesUnknownGroupReferenceParser1)
{
    ASSERT_NO_THROW_GMX(sc_.setIndexGroups(nullptr));
    EXPECT_THROW_GMX(sc_.parseFromString("group \"foo\""), gmx::InconsistentInputError);
    EXPECT_THROW_GMX(sc_.parseFromString("4"), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, HandlesUnknownGroupReferenceParser2)
{
    ASSERT_NO_THROW_GMX(loadIndexGroups("simple.ndx"));
    EXPECT_THROW_GMX(sc_.parseFromString("group \"foo\""), gmx::InconsistentInputError);
    EXPECT_THROW_GMX(sc_.parseFromString("4"), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, HandlesUnknownGroupReferenceDelayed1)
{
    ASSERT_NO_THROW_GMX(sc_.parseFromString("group \"foo\""));
    ASSERT_NO_FATAL_FAILURE(setAtomCount(10));
    EXPECT_THROW_GMX(sc_.setIndexGroups(nullptr), gmx::InconsistentInputError);
    EXPECT_THROW_GMX(sc_.compile(), gmx::APIError);
}

TEST_F(SelectionCollectionTest, HandlesUnknownGroupReferenceDelayed2)
{
    ASSERT_NO_THROW_GMX(sc_.parseFromString("group 4; group \"foo\""));
    ASSERT_NO_FATAL_FAILURE(setAtomCount(10));
    EXPECT_THROW_GMX(loadIndexGroups("simple.ndx"), gmx::InconsistentInputError);
    EXPECT_THROW_GMX(sc_.compile(), gmx::APIError);
}

TEST_F(SelectionCollectionTest, HandlesUnsortedGroupReference)
{
    ASSERT_NO_THROW_GMX(loadIndexGroups("simple.ndx"));
    EXPECT_THROW_GMX(sc_.parseFromString("atomnr 1 to 3 and group \"GrpUnsorted\""),
                     gmx::InconsistentInputError);
    EXPECT_THROW_GMX(sc_.parseFromString("group 2 or atomnr 2 to 5"),
                     gmx::InconsistentInputError);
    EXPECT_THROW_GMX(sc_.parseFromString("within 1 of group 2"),
                     gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, HandlesUnsortedGroupReferenceDelayed)
{
    ASSERT_NO_THROW_GMX(sc_.parseFromString("atomnr 1 to 3 and group \"GrpUnsorted\""));
    ASSERT_NO_THROW_GMX(sc_.parseFromString("atomnr 1 to 3 and group 2"));
    EXPECT_THROW_GMX(loadIndexGroups("simple.ndx"), gmx::InconsistentInputError);
    // TODO: Add a separate check in the selection compiler for a safer API
    // (makes sense in the future if the compiler needs the information for
    // other purposes as well).
    // EXPECT_THROW_GMX(sc_.compile(), gmx::APIError);
}

TEST_F(SelectionCollectionTest, HandlesOutOfRangeAtomIndexInGroup)
{
    ASSERT_NO_THROW_GMX(sc_.setTopology(nullptr, 5));
    ASSERT_NO_THROW_GMX(loadIndexGroups("simple.ndx"));
    EXPECT_THROW_GMX(sc_.parseFromString("group \"GrpB\""), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, HandlesOutOfRangeAtomIndexInGroupDelayed)
{
    ASSERT_NO_THROW_GMX(loadIndexGroups("simple.ndx"));
    ASSERT_NO_THROW_GMX(sc_.parseFromString("group \"GrpB\""));
    EXPECT_THROW_GMX(sc_.setTopology(nullptr, 5), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, HandlesOutOfRangeAtomIndexInGroupDelayed2)
{
    ASSERT_NO_THROW_GMX(sc_.setTopology(nullptr, 5));
    ASSERT_NO_THROW_GMX(sc_.parseFromString("group \"GrpB\""));
    EXPECT_THROW_GMX(loadIndexGroups("simple.ndx"), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromMissingMoleculeInfo)
{
    ASSERT_NO_THROW_GMX(sc_.parseFromString("molindex 1 to 5"));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW_GMX(sc_.compile(), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromMissingAtomTypes)
{
    ASSERT_NO_THROW_GMX(sc_.parseFromString("type CA"));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW_GMX(sc_.compile(), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromMissingPDBInfo)
{
    ASSERT_NO_THROW_GMX(sc_.parseFromString("altloc A"));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    EXPECT_THROW_GMX(sc_.compile(), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromInvalidPermutation)
{
    ASSERT_NO_THROW_GMX(sc_.parseFromString("all permute 1 1"));
    ASSERT_NO_FATAL_FAILURE(setAtomCount(10));
    EXPECT_THROW_GMX(sc_.compile(), gmx::InvalidInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromInvalidPermutation2)
{
    ASSERT_NO_THROW_GMX(sc_.parseFromString("all permute 3 2 1"));
    ASSERT_NO_FATAL_FAILURE(setAtomCount(10));
    EXPECT_THROW_GMX(sc_.compile(), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, RecoversFromInvalidPermutation3)
{
    ASSERT_NO_THROW_GMX(sc_.parseFromString("x < 1.5 permute 3 2 1"));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    ASSERT_NO_THROW_GMX(sc_.compile());
    EXPECT_THROW_GMX(sc_.evaluate(topManager_.frame(), nullptr), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, HandlesFramesWithTooSmallAtomSubsets)
{
    ASSERT_NO_THROW_GMX(sc_.parseFromString("atomnr 3 to 10"));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    ASSERT_NO_THROW_GMX(sc_.compile());
    topManager_.frame()->natoms = 8;
    EXPECT_THROW_GMX(sc_.evaluate(topManager_.frame(), nullptr), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, HandlesFramesWithTooSmallAtomSubsets2)
{
    const int index[] = { 1, 2, 3, 9 };
    ASSERT_NO_THROW_GMX(sc_.parseFromString("atomnr 3 4 7 10"));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    ASSERT_NO_THROW_GMX(sc_.compile());
    topManager_.initFrameIndices(index);
    EXPECT_THROW_GMX(sc_.evaluate(topManager_.frame(), nullptr), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, HandlesFramesWithTooSmallAtomSubsets3)
{
    const int index[] = { 0, 1, 2, 3, 4, 5, 6, 9, 10, 11 };
    // Evaluating the positions will require atoms 1-3, 7-12.
    ASSERT_NO_THROW_GMX(sc_.parseFromString("whole_res_cog of atomnr 2 7 11"));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    ASSERT_NO_THROW_GMX(sc_.compile());
    topManager_.initFrameIndices(index);
    EXPECT_THROW_GMX(sc_.evaluate(topManager_.frame(), nullptr), gmx::InconsistentInputError);
}

TEST_F(SelectionCollectionTest, HandlesFramesWithTooSmallAtomSubsets4)
{
    ASSERT_NO_THROW_GMX(sc_.parseFromString("mindistance from atomnr 1 to 5 < 2"));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    ASSERT_NO_THROW_GMX(sc_.compile());
    topManager_.frame()->natoms = 10;
    EXPECT_THROW_GMX(sc_.evaluate(topManager_.frame(), nullptr), gmx::InconsistentInputError);
}

// TODO: Tests for more evaluation errors

/********************************************************************
 * Tests for interactive selection input
 */

TEST_F(SelectionCollectionInteractiveTest, HandlesBasicInput)
{
    const char *const input[] = {
        "foo = resname RA",
        "resname RB",
        "\"Name\" resname RC"
    };
    runTest(-1, true, input);
}

TEST_F(SelectionCollectionInteractiveTest, HandlesContinuation)
{
    const char *const input[] = {
        "resname RB and \\",
        "resname RC"
    };
    runTest(-1, true, input);
}

TEST_F(SelectionCollectionInteractiveTest, HandlesSingleSelectionInput)
{
    const char *const input[] = {
        "foo = resname RA",
        "resname RA"
    };
    runTest(1, true, input);
}

TEST_F(SelectionCollectionInteractiveTest, HandlesTwoSelectionInput)
{
    const char *const input[] = {
        "resname RA",
        "resname RB"
    };
    runTest(2, true, input);
}

TEST_F(SelectionCollectionInteractiveTest, HandlesStatusWithGroups)
{
    const char *const input[] = {
        "resname RA",
        ""
    };
    loadIndexGroups("simple.ndx");
    runTest(-1, true, input);
}

TEST_F(SelectionCollectionInteractiveTest, HandlesStatusWithExistingSelections)
{
    const char *const input[] = {
        "",
        "bar = resname RC",
        "resname RA",
        ""
    };
    ASSERT_NO_THROW_GMX(sc_.parseFromString("foo = resname RA"));
    ASSERT_NO_THROW_GMX(sc_.parseFromString("resname RB"));
    runTest(-1, true, input);
}

TEST_F(SelectionCollectionInteractiveTest, HandlesSingleSelectionInputStatus)
{
    const char *const input[] = {
        "foo = resname RA",
        "",
        "resname RB"
    };
    runTest(1, true, input);
}

TEST_F(SelectionCollectionInteractiveTest, HandlesTwoSelectionInputStatus)
{
    const char *const input[] = {
        "\"Sel\" resname RA",
        "",
        "resname RB"
    };
    runTest(2, true, input);
}

TEST_F(SelectionCollectionInteractiveTest, HandlesMultiSelectionInputStatus)
{
    const char *const input[] = {
        "\"Sel\" resname RA",
        "\"Sel2\" resname RB",
        ""
    };
    runTest(-1, true, input);
}

TEST_F(SelectionCollectionInteractiveTest, HandlesNoFinalNewline)
{
    // TODO: There is an extra prompt printed after the input is finished; it
    // would be cleaner not to have it, but it's only a cosmetic issue.
    const char *const input[] = {
        "resname RA"
    };
    helper_.setLastNewline(false);
    runTest(-1, true, input);
}

TEST_F(SelectionCollectionInteractiveTest, HandlesEmptySelections)
{
    const char *const input[] = {
        "resname RA;",
        "; resname RB;;",
        " ",
        ";"
    };
    runTest(-1, true, input);
}

TEST_F(SelectionCollectionInteractiveTest, HandlesMultipleSelectionsOnLine)
{
    const char *const input[] = {
        "resname RA; resname RB and \\",
        "resname RC"
    };
    runTest(2, true, input);
}

TEST_F(SelectionCollectionInteractiveTest, HandlesNoninteractiveInput)
{
    const char *const input[] = {
        "foo = resname RA",
        "resname RB",
        "\"Name\" resname RC"
    };
    runTest(-1, false, input);
}

TEST_F(SelectionCollectionInteractiveTest, HandlesSingleSelectionInputNoninteractively)
{
    const char *const input[] = {
        "foo = resname RA",
        "resname RA"
    };
    runTest(1, false, input);
}


/********************************************************************
 * Tests for selection keywords
 */

TEST_F(SelectionCollectionDataTest, HandlesAllNone)
{
    static const char * const selections[] = {
        "all",
        "none"
    };
    runTest(10, selections);
}

TEST_F(SelectionCollectionDataTest, HandlesAtomnr)
{
    static const char * const selections[] = {
        "atomnr 1 to 3 6 to 8",
        "atomnr 4 2 5 to 7",
        "atomnr <= 5"
    };
    runTest(10, selections);
}

TEST_F(SelectionCollectionDataTest, HandlesResnr)
{
    static const char * const selections[] = {
        "resnr 1 2 5",
        "resid 4 to 3"
    };
    runTest("simple.gro", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesResIndex)
{
    static const char * const selections[] = {
        "resindex 1 4",
        "residue 1 3"
    };
    runTest("simple.pdb", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesMolIndex)
{
    static const char * const selections[] = {
        "molindex 1 4",
        "molecule 2 3 5"
    };
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    ASSERT_NO_FATAL_FAILURE(topManager_.loadTopology("simple.gro"));
    topManager_.initUniformMolecules(3);
    ASSERT_NO_FATAL_FAILURE(setTopology());
    ASSERT_NO_FATAL_FAILURE(runCompiler());
}

TEST_F(SelectionCollectionDataTest, HandlesAtomname)
{
    static const char * const selections[] = {
        "name CB",
        "atomname S1 S2"
    };
    runTest("simple.gro", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesPdbAtomname)
{
    static const char * const selections[] = {
        "name HG21",
        "name 1HG2",
        "pdbname HG21 CB",
        "pdbatomname 1HG2"
    };
    runTest("simple.pdb", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesAtomtype)
{
    static const char * const selections[] = {
        "atomtype CA"
    };
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    ASSERT_NO_FATAL_FAILURE(topManager_.loadTopology("simple.gro"));
    const char *const types[] = { "CA", "SA", "SB" };
    topManager_.initAtomTypes(types);
    ASSERT_NO_FATAL_FAILURE(setTopology());
    ASSERT_NO_FATAL_FAILURE(runCompiler());
}

TEST_F(SelectionCollectionDataTest, HandlesChain)
{
    static const char * const selections[] = {
        "chain A",
        "chain B"
    };
    runTest("simple.pdb", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesMass)
{
    static const char * const selections[] = {
        "mass > 5"
    };
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    EXPECT_TRUE(sc_.requiredTopologyProperties().needsMasses);
    ASSERT_NO_FATAL_FAILURE(topManager_.loadTopology("simple.gro"));
    t_atoms &atoms = topManager_.atoms();
    for (int i = 0; i < atoms.nr; ++i)
    {
        atoms.atom[i].m = 1.0 + i;
    }
    atoms.haveMass = TRUE;
    ASSERT_NO_FATAL_FAILURE(setTopology());
    ASSERT_NO_FATAL_FAILURE(runCompiler());
}

TEST_F(SelectionCollectionDataTest, HandlesCharge)
{
    static const char * const selections[] = {
        "charge < 0.5"
    };
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    ASSERT_NO_FATAL_FAILURE(topManager_.loadTopology("simple.gro"));
    t_atoms &atoms = topManager_.atoms();
    for (int i = 0; i < atoms.nr; ++i)
    {
        atoms.atom[i].q = i / 10.0;
    }
    // Ensure exact representation of 0.5 is used, so that the test is
    // reproducible.
    atoms.atom[5].q  = 0.5;
    atoms.haveCharge = TRUE;
    ASSERT_NO_FATAL_FAILURE(setTopology());
    ASSERT_NO_FATAL_FAILURE(runCompiler());
}

TEST_F(SelectionCollectionDataTest, HandlesAltLoc)
{
    static const char * const selections[] = {
        "altloc \" \"",
        "altloc A"
    };
    runTest("simple.pdb", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesInsertCode)
{
    static const char * const selections[] = {
        "insertcode \" \"",
        "insertcode A"
    };
    runTest("simple.pdb", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesOccupancy)
{
    static const char * const selections[] = {
        "occupancy 1",
        "occupancy < .5"
    };
    runTest("simple.pdb", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesBeta)
{
    static const char * const selections[] = {
        "beta 0",
        "beta >= 0.3"
    };
    runTest("simple.pdb", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesResname)
{
    static const char * const selections[] = {
        "resname RA",
        "resname RB RC"
    };
    runTest("simple.gro", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesCoordinateKeywords)
{
    static const char * const selections[] = {
        "x < 3",
        "y >= 3",
        "x {-1 to 2}"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesSameResidue)
{
    static const char * const selections[] = {
        "same residue as atomnr 1 4 12"
    };
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesSameResidueName)
{
    static const char * const selections[] = {
        "same resname as atomnr 1 14"
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
        "dyn_res_cog of x < 3"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates
             | efTestPositionAtoms);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesDistanceKeyword)
{
    static const char * const selections[] = {
        "distance from cog of resnr 1 < 2"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesMinDistanceKeyword)
{
    static const char * const selections[] = {
        "mindistance from resnr 1 < 2"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesWithinKeyword)
{
    static const char * const selections[] = {
        "within 1 of resnr 2"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesInSolidAngleKeyword)
{
    // Both of these should evaluate to empty on a correct implementation.
    static const char * const selections[] = {
        "resname TP and not insolidangle center cog of resname C span resname R cutoff 20",
        "resname TN and insolidangle center cog of resname C span resname R cutoff 20"
    };
    setFlags(TestFlags() | efDontTestCompiledAtoms | efTestEvaluation);
    runTest("sphere.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesPermuteModifier)
{
    static const char * const selections[] = {
        "all permute 3 1 2",
        "res_cog of resnr 1 to 4 permute 2 1",
        "name CB S1 and res_cog x < 3 permute 2 1"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates
             | efTestPositionAtoms | efTestPositionMapping);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesPlusModifier)
{
    static const char * const selections[] = {
        "name S2 plus name S1",
        "res_cog of resnr 2 plus res_cog of resnr 1 plus res_cog of resnr 3",
        "name S1 and y < 3 plus res_cog of x < 2.5"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates
             | efTestPositionAtoms | efTestPositionMapping);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesMergeModifier)
{
    static const char * const selections[] = {
        "name S2 merge name S1",
        "resnr 1 2 and name S2 merge resnr 1 2 and name S1 merge res_cog of resnr 1 2",
        "name S1 and x < 2.5 merge res_cog of x < 2.5"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates
             | efTestPositionAtoms | efTestPositionMapping);
    runTest("simple.gro", selections);
}


/********************************************************************
 * Tests for generic selection evaluation
 */

TEST_F(SelectionCollectionDataTest, ComputesMassesAndCharges)
{
    static const char * const selections[] = {
        "name CB",
        "y > 2",
        "res_cog of y > 2"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionAtoms
             | efTestPositionMasses | efTestPositionCharges);
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    ASSERT_NO_FATAL_FAILURE(topManager_.loadTopology("simple.gro"));
    t_atoms &atoms = topManager_.atoms();
    for (int i = 0; i < atoms.nr; ++i)
    {
        atoms.atom[i].m =   1.0 + i / 100.0;
        atoms.atom[i].q = -(1.0 + i / 100.0);
    }
    atoms.haveMass   = TRUE;
    atoms.haveCharge = TRUE;
    ASSERT_NO_FATAL_FAILURE(setTopology());
    ASSERT_NO_FATAL_FAILURE(runCompiler());
    ASSERT_NO_FATAL_FAILURE(runEvaluate());
    ASSERT_NO_FATAL_FAILURE(runEvaluateFinal());
}

TEST_F(SelectionCollectionDataTest, ComputesMassesAndChargesWithoutTopology)
{
    static const char * const selections[] = {
        "atomnr 1 to 3 8 to 9",
        "y > 2",
        "cog of (y > 2)"
    };
    setFlags(TestFlags() | efTestPositionAtoms
             | efTestPositionMasses | efTestPositionCharges);
    runTest(10, selections);
}

TEST_F(SelectionCollectionDataTest, HandlesFramesWithAtomSubsets)
{
    const int          index[]      = { 0, 1, 2, 3, 4, 5, 9, 10, 11 };
    const char * const selections[] = {
        "resnr 1 4",
        "atomnr 1 2 5 11 and y > 2",
        "res_cog of atomnr 2 5 11"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionAtoms);
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    ASSERT_NO_FATAL_FAILURE(runCompiler());
    topManager_.initFrameIndices(index);
    ASSERT_NO_FATAL_FAILURE(runEvaluate());
    ASSERT_NO_FATAL_FAILURE(runEvaluateFinal());
}


/********************************************************************
 * Tests for selection syntactic constructs
 */

TEST_F(SelectionCollectionDataTest, HandlesSelectionNames)
{
    static const char * const selections[] = {
        "\"GroupSelection\" group \"GrpA\"",
        "\"DynamicSelection\" x < 5",
        "y < 3"
    };
    setFlags(TestFlags() | efTestSelectionNames);
    ASSERT_NO_THROW_GMX(loadIndexGroups("simple.ndx"));
    runTest(10, selections);
}

TEST_F(SelectionCollectionDataTest, HandlesIndexGroupsInSelections)
{
    static const char * const selections[] = {
        "group \"GrpA\"",
        "GrpB",
        "1",
        // These test that the name of the group is not too eagerly promoted to
        // the name of the selection.
        "group \"GrpB\" and resname RB",
        "group \"GrpA\" permute 5 3 2 1 4",
        "group \"GrpA\" plus group \"GrpB\"",
        "res_cog of group \"GrpA\""
    };
    setFlags(TestFlags() | efTestSelectionNames);
    ASSERT_NO_THROW_GMX(loadIndexGroups("simple.ndx"));
    runTest("simple.gro", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesIndexGroupsInSelectionsDelayed)
{
    static const char * const selections[] = {
        "group \"GrpA\"",
        "GrpB",
        "1",
        "group \"GrpB\" and resname RB"
    };
    setFlags(TestFlags() | efTestSelectionNames);
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    ASSERT_NO_THROW_GMX(loadIndexGroups("simple.ndx"));
    ASSERT_NO_FATAL_FAILURE(runCompiler());
}

TEST_F(SelectionCollectionDataTest, HandlesUnsortedIndexGroupsInSelections)
{
    static const char * const selections[] = {
        "foo = group \"GrpUnsorted\"",
        "group \"GrpUnsorted\"",
        "GrpUnsorted",
        "2",
        "res_cog of group \"GrpUnsorted\"",
        "group \"GrpUnsorted\" permute 2 1",
        "foo"
    };
    setFlags(TestFlags() | efTestPositionAtoms | efTestPositionMapping
             | efTestSelectionNames);
    ASSERT_NO_THROW_GMX(loadIndexGroups("simple.ndx"));
    runTest("simple.gro", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesUnsortedIndexGroupsInSelectionsDelayed)
{
    static const char * const selections[] = {
        "foo = group \"GrpUnsorted\"",
        "group \"GrpUnsorted\"",
        "GrpUnsorted",
        "2",
        "res_cog of group \"GrpUnsorted\"",
        "group \"GrpUnsorted\" permute 2 1",
        "foo"
    };
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    ASSERT_NO_FATAL_FAILURE(loadTopology("simple.gro"));
    ASSERT_NO_THROW_GMX(loadIndexGroups("simple.ndx"));
    ASSERT_NO_FATAL_FAILURE(runCompiler());
}


TEST_F(SelectionCollectionDataTest, HandlesConstantPositions)
{
    static const char * const selections[] = {
        "[1, -2, 3.5]"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates
             | efTestPositionMapping);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesConstantPositionsWithModifiers)
{
    static const char * const selections[] = {
        "[0, 0, 0] plus [0, 1, 0]"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates
             | efTestPositionMapping);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesWithinConstantPositions)
{
    static const char * const selections[] = {
        "within 1 of [2, 1, 0]"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesOverlappingIntegerRanges)
{
    static const char * const selections[] = {
        "atomnr 2 to 4 5 to 8",
        "atomnr 2 to 5 4 to 7"
    };
    ASSERT_NO_FATAL_FAILURE(runTest(10, selections));
}


TEST_F(SelectionCollectionDataTest, HandlesOverlappingRealRanges)
{
    static const char * const selections[] = {
        "charge {-0.35 to -0.05 0.25 to 0.75}",
        "charge {0.05 to -0.3 -0.05 to 0.55}"
    };
    ASSERT_NO_FATAL_FAILURE(runParser(selections));
    ASSERT_NO_FATAL_FAILURE(topManager_.loadTopology("simple.gro"));
    t_atoms &atoms = topManager_.atoms();
    for (int i = 0; i < atoms.nr; ++i)
    {
        atoms.atom[i].q = i / 10.0 - 0.5;
    }
    atoms.haveCharge = TRUE;
    ASSERT_NO_FATAL_FAILURE(setTopology());
    ASSERT_NO_FATAL_FAILURE(runCompiler());
}


TEST_F(SelectionCollectionDataTest, HandlesForcedStringMatchingMode)
{
    static const char * const selections[] = {
        "name = S1 \"C?\"",
        "name ? S1 \"C?\""
    };
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesWildcardMatching)
{
    static const char * const selections[] = {
        "name \"S?\"",
        "name ? \"S?\""
    };
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesRegexMatching)
{
    static const char * const selections[] = {
        "resname \"R[BD]\"",
        "resname ~ \"R[BD]\""
    };
    if (gmx::Regex::isSupported())
    {
        runTest("simple.gro", selections);
    }
}


TEST_F(SelectionCollectionDataTest, HandlesBasicBoolean)
{
    static const char * const selections[] = {
        "atomnr 1 to 5 and atomnr 2 to 7",
        "atomnr 1 to 5 or not atomnr 3 to 8",
        "not not atomnr 1 to 5 and atomnr 2 to 6 and not not atomnr 3 to 7",
        "atomnr 1 to 5 and (atomnr 2 to 7 and atomnr 3 to 6)",
        "x < 5 and atomnr 1 to 5 and y < 3 and atomnr 2 to 4"
    };
    runTest(10, selections);
}


TEST_F(SelectionCollectionDataTest, HandlesDynamicAtomValuedParameters)
{
    static const char * const selections[] = {
        "same residue as (atomnr 3 5 13 or y > 5)",
        "(resnr 1 3 5 or x > 10) and same residue as (atomnr 3 5 13 or z > 5)"
    };
    setFlags(TestFlags() | efTestEvaluation);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesEmptySelectionWithUnevaluatedExpressions)
{
    static const char * const selections[] = {
        "none and x > 2",
        "none and same resname as resnr 2"
    };
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesEmptyReferenceForSame)
{
    static const char * const selections[] = {
        "same residue as none",
        "same resname as none"
    };
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesPositionModifiersForKeywords)
{
    static const char * const selections[] = {
        "res_cog x > 2",
        "name CB and res_cog y > 2.5"
    };
    setFlags(TestFlags() | efTestEvaluation);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesPositionModifiersForMethods)
{
    static const char * const selections[] = {
        "res_cog distance from cog of resnr 1 < 2",
        "res_cog within 2 of cog of resnr 1"
    };
    setFlags(TestFlags() | efTestEvaluation);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesKeywordOfPositions)
{
    static const char * const selections[] = {
        "x < y of cog of resnr 2"
    };
    setFlags(TestFlags() | efTestEvaluation);
    runTest("simple.gro", selections);
}

TEST_F(SelectionCollectionDataTest, HandlesKeywordOfPositionsInArithmetic)
{
    static const char * const selections[] = {
        "x - y of cog of resnr 2 < 0"
    };
    setFlags(TestFlags() | efTestEvaluation);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesNumericComparisons)
{
    static const char * const selections[] = {
        "x > 2",
        "2 < x",
        "y > resnr",
        "resnr < 2.5",
        "2.5 > resnr"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesArithmeticExpressions)
{
    static const char * const selections[] = {
        "x+1 > 3",
        "(y-1)^2 <= 1",
        "x+--1 > 3",
        "-x+-1 < -3"
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
        "index < 3"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesComplexNumericVariables)
{
    static const char * const selections[] = {
        "value = x + y",
        "resname RA and value <= 4",
        "resname RA RB and x < 3 and value <= 4",
        "index = atomnr",
        "resname RA and index < 3",
        "resname RB and y < 3 and index < 6"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesPositionVariables)
{
    static const char * const selections[] = {
        "foo = res_cog of resname RA",
        "foo",
        "within 1 of foo",
        "bar = cog of resname RA",
        "bar",
        "within 1 of bar"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesPositionVariableInModifier)
{
    static const char * const selections[] = {
        "foo = cog of resnr 1",
        "cog of resnr 2 plus foo",
        "cog of resnr 3 plus foo"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesConstantPositionInVariable)
{
    static const char * const selections[] = {
        "constpos = [1.0, 2.5, 0.5]",
        "constpos",
        "within 2 of constpos"
    };
    setFlags(TestFlags() | efTestEvaluation | efTestPositionCoordinates
             | efTestPositionAtoms);
    runTest("simple.gro", selections);
}


TEST_F(SelectionCollectionDataTest, HandlesNumericConstantsInVariables)
{
    static const char * const selections[] = {
        "constint = 4",
        "constreal1 = 0.5",
        "constreal2 = 2.7",
        "resnr < constint",
        "x + constreal1 < constreal2"
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
        "atomnr 1 to 5 or (atomnr 4 to 6 and (atomnr 5 to 7 or x < 2))"
    };
    runTest(10, selections);
}


TEST_F(SelectionCollectionDataTest, HandlesBooleanStaticAnalysisWithVariables)
{
    static const char * const selections[] = {
        "foo = atomnr 4 to 7 or x < 2",
        "atomnr 1 to 4 and foo",
        "atomnr 2 to 6 and y < 3 and foo",
        "atomnr 6 to 10 and not foo"
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
        "atomnr 6 to 10 and not foo"
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
        "unused2 = atomnr 3 to 5"
    };
    runTest(10, selections);
}


TEST_F(SelectionCollectionDataTest, HandlesVariablesWithStaticEvaluationGroups)
{
    static const char * const selections[] = {
        "foo = atomnr 4 to 7 and x < 2",
        "atomnr 1 to 5 and foo",
        "atomnr 3 to 7 and foo"
    };
    runTest(10, selections);
}


TEST_F(SelectionCollectionDataTest, HandlesVariablesWithMixedEvaluationGroups)
{
    static const char * const selections[] = {
        "foo = atomnr 4 to 7 and x < 2",
        "atomnr 1 to 6 and foo",
        "within 1 of foo",
        "foo"
    };
    runTest(10, selections);
}


TEST_F(SelectionCollectionDataTest, HandlesVariablesWithMixedEvaluationGroups2)
{
    static const char * const selections[] = {
        "foo = atomnr 1 to 8 and x < 10",
        "atomnr 1 to 5 and y < 10 and foo",
        "foo"
    };
    setFlags(TestFlags() | efTestEvaluation);
    runTest("simple.gro", selections);
}


} // namespace
