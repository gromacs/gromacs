/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * \brief Implemention of functions for comparing trajectories
 * produced by mdrun.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "trajectorycomparison.h"

#include <array>
#include <iterator>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/strconvert.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"
#include "testutils/trajectoryreader.h"

#include "programs/mdrun/tests/comparison_helpers.h"

namespace gmx
{
namespace test
{

using ::testing::Pointwise;

/*! \brief Compares the box from \c reference and \c test
 * according to the \c matchSettings and \c tolerance.
 *
 * \todo This could be streamlined when we have a proper 3D matrix
 * class and view. */
static void compareBox(const TrajectoryFrame&              reference,
                       const TrajectoryFrame&              test,
                       const TrajectoryFrameMatchSettings& matchSettings,
                       const FloatingPointTolerance        tolerance)
{
    if (!matchSettings.mustCompareBox)
    {
        return;
    }
    bool canCompareBox = true;
    if (!reference.hasBox())
    {
        ADD_FAILURE() << "Comparing the box was required, "
                         "but the reference frame did not have one";
        canCompareBox = false;
    }
    if (!test.hasBox())
    {
        ADD_FAILURE() << "Comparing the box was required, "
                         "but the test frame did not have one";
        canCompareBox = false;
    }
    if (!canCompareBox)
    {
        return;
    }

    // Do the comparing.
    for (int d = 0; d < DIM; ++d)
    {
        for (int dd = 0; dd < DIM; ++dd)
        {
            EXPECT_REAL_EQ_TOL(reference.box()[d][dd], test.box()[d][dd], tolerance);
        }
    }
}

/*! \brief Help put all atom coordinates in \c frame into its box.
 *
 * This can perhaps go away when frame->x is a container. */
static std::vector<RVec> putAtomsInBox(const TrajectoryFrame& frame)
{
    std::vector<RVec> x(frame.x().begin(), frame.x().end());
    matrix            box;
    for (int d = 0; d < DIM; ++d)
    {
        for (int dd = 0; dd < DIM; ++dd)
        {
            box[d][dd] = frame.box()[d][dd];
        }
    }
    // Note we don't need to compare bPBC because put_atoms_in_box
    // implements a fallback if nothing specific was set in the
    // trajectory frame.
    put_atoms_in_box(frame.pbc(), box, x);
    return x;
}

/*! \brief Return whether the \c comparisonConditions and emptiness of
 * reference and test frames means that a comparison should be
 * attempted.
 *
 * This allows the framework to determine whether it is an error if a
 * comparison cannot be made. */
static bool shouldDoComparison(const ComparisonConditions comparisonConditions,
                               const bool                 referenceIsEmpty,
                               const bool                 testIsEmpty)
{
    if (comparisonConditions == ComparisonConditions::NoComparison)
    {
        return false;
    }

    bool doComparison = true;
    if (testIsEmpty)
    {
        if (comparisonConditions == ComparisonConditions::MustCompare
            || comparisonConditions == ComparisonConditions::CompareIfBothFound
            || (!referenceIsEmpty && comparisonConditions == ComparisonConditions::CompareIfReferenceFound))
        {
            ADD_FAILURE() << "Test frame lacked quantity for required comparison";
        }
        doComparison = false;
    }
    if (referenceIsEmpty)
    {
        if (comparisonConditions == ComparisonConditions::MustCompare
            || comparisonConditions == ComparisonConditions::CompareIfBothFound
            || (!testIsEmpty && comparisonConditions == ComparisonConditions::CompareIfTestFound))
        {
            ADD_FAILURE() << "Reference frame lacked quantity for required comparison";
        }
        doComparison = false;
    }
    return doComparison;
}

/*! \brief Compares the position coordinates from \c reference and \c test
 * according to the \c matchSettings and \c tolerance. */
static void compareCoordinates(const TrajectoryFrame&              reference,
                               const TrajectoryFrame&              test,
                               const TrajectoryFrameMatchSettings& matchSettings,
                               const FloatingPointTolerance        tolerance)
{
    SCOPED_TRACE("Comparing coordinates");
    if (!shouldDoComparison(matchSettings.coordinatesComparison, reference.x().empty(), test.x().empty()))
    {
        return;
    }

    bool canHandlePbc = true;
    if (!reference.hasBox())
    {
        ADD_FAILURE() << "Comparing positions required PBC handling, "
                         "but the reference frame did not have a box";
        canHandlePbc = false;
    }
    if (!test.hasBox())
    {
        ADD_FAILURE() << "Comparing positions required PBC handling, "
                         "but the test frame did not have a box";
        canHandlePbc = false;
    }

    if (matchSettings.requirePbcHandling && !canHandlePbc)
    {
        ADD_FAILURE() << "Cannot compare positions for the above reason(s)";
        return;
    }

    if ((matchSettings.handlePbcIfPossible || matchSettings.requirePbcHandling) && canHandlePbc)
    {
        EXPECT_THAT(putAtomsInBox(test), Pointwise(RVecEq(tolerance), putAtomsInBox(reference)));
    }
    else
    {
        EXPECT_THAT(test.x(), Pointwise(RVecEq(tolerance), reference.x()));
    }
}

/*! \brief Compares the velocities from \c reference and \c test
 * according to the \c matchSettings and \c tolerance. */
static void compareVelocities(const TrajectoryFrame&              reference,
                              const TrajectoryFrame&              test,
                              const TrajectoryFrameMatchSettings& matchSettings,
                              const FloatingPointTolerance        tolerance)
{
    SCOPED_TRACE("Comparing velocities");
    if (!shouldDoComparison(matchSettings.velocitiesComparison, reference.v().empty(), test.v().empty()))
    {
        return;
    }
    EXPECT_THAT(test.v(), Pointwise(RVecEq(tolerance), reference.v()));
}

/*! \brief Compares the forces from \c reference and \c test
 * according to the \c matchSettings and \c tolerance. */
static void compareForces(const TrajectoryFrame&              reference,
                          const TrajectoryFrame&              test,
                          const TrajectoryFrameMatchSettings& matchSettings,
                          const FloatingPointTolerance        tolerance)
{
    SCOPED_TRACE("Comparing forces");
    if (!shouldDoComparison(matchSettings.forcesComparison, reference.f().empty(), test.f().empty()))
    {
        return;
    }
    EXPECT_THAT(test.f(), Pointwise(RVecEq(tolerance), reference.f()));
}


const TrajectoryTolerances TrajectoryComparison::s_defaultTrajectoryTolerances{
    defaultRealTolerance(),                                               // box
    relativeToleranceAsFloatingPoint(1.0, 1.0e-3),                        // positions
    defaultRealTolerance(),                                               // velocities
    relativeToleranceAsFloatingPoint(100.0, GMX_DOUBLE ? 1.0e-7 : 5.0e-5) // forces
};


TrajectoryComparison::TrajectoryComparison(const TrajectoryFrameMatchSettings& matchSettings,
                                           const TrajectoryTolerances&         tolerances) :
    matchSettings_(matchSettings), tolerances_(tolerances)
{
}

void TrajectoryComparison::operator()(const TrajectoryFrame& reference, const TrajectoryFrame& test) const
{
    if (numComparedFrames_ >= matchSettings_.maxNumTrajectoryFrames)
    {
        // Nothing should be compared
        return;
    }
    SCOPED_TRACE("Comparing trajectory reference frame " + reference.frameName()
                 + " and test frame " + test.frameName());
    EXPECT_EQ(reference.step(), test.step());
    EXPECT_EQ(reference.time(), test.time());
    compareBox(reference, test, matchSettings_, tolerances_.box);
    compareCoordinates(reference, test, matchSettings_, tolerances_.coordinates);
    compareVelocities(reference, test, matchSettings_, tolerances_.velocities);
    compareForces(reference, test, matchSettings_, tolerances_.forces);
    numComparedFrames_++;
}

/*! \brief Return whether the \c comparisonConditions and emptiness of
 * the test frame means that a comparison should be attempted.
 *
 * This allows the framework to determine whether it is an error if a
 * comparison cannot be made. */
static bool shouldDoComparison(ArrayRef<const RVec> values, ComparisonConditions comparisonConditions)
{
    if (comparisonConditions == ComparisonConditions::NoComparison)
    {
        return false;
    }
    if (values.empty())
    {
        // Empty array ref is only acceptable if comparison is not required
        if (comparisonConditions != ComparisonConditions::CompareIfTestFound
            && comparisonConditions != ComparisonConditions::CompareIfBothFound)
        {
            ADD_FAILURE() << "Test frame lacked quantity for required comparison";
        }
        return false;
    }
    return true;
}

/*! \brief Compares the positions from \c frame to reference data
 * according to the \c matchSettings and \c tolerance.
 */
static void checkPositionsAgainstReference(const TrajectoryFrame&              frame,
                                           const TrajectoryFrameMatchSettings& matchSettings,
                                           const TrajectoryTolerances&         tolerance,
                                           TestReferenceChecker*               checker)
{
    SCOPED_TRACE("Comparing positions");
    if (shouldDoComparison(frame.x(), matchSettings.coordinatesComparison))
    {
        ArrayRef<const RVec> positions{};
        std::vector<RVec>    shiftedPositions{};
        if (frame.hasBox() && (matchSettings.handlePbcIfPossible || matchSettings.requirePbcHandling))
        {
            shiftedPositions = putAtomsInBox(frame);
            positions        = shiftedPositions;
        }
        else
        {
            positions = frame.x();
        }

        if (!frame.hasBox() && matchSettings.requirePbcHandling)
        {
            ADD_FAILURE() << "Comparing positions required PBC handling, "
                             "but the test frame did not have a box";
            return;
        }
        checker->setDefaultTolerance(tolerance.coordinates);
        std::string id = frame.frameName() + " X";
        checker->checkSequence(std::begin(positions), std::end(positions), id.c_str());
    }
}

/*! \brief Compares the velocities from \c frame to reference data
 * according to the \c matchSettings and \c tolerance.
 */
static void checkVelocitiesAgainstReference(const TrajectoryFrame&              frame,
                                            const TrajectoryFrameMatchSettings& matchSettings,
                                            const TrajectoryTolerances&         tolerance,
                                            TestReferenceChecker*               checker)
{
    SCOPED_TRACE("Comparing velocities");
    if (shouldDoComparison(frame.v(), matchSettings.velocitiesComparison))
    {
        checker->setDefaultTolerance(tolerance.velocities);
        std::string id = frame.frameName() + " V";
        checker->checkSequence(std::begin(frame.v()), std::end(frame.v()), id.c_str());
    }
}

/*! \brief Compares the forces from \c frame to reference data
 * according to the \c matchSettings and \c tolerance.
 */
static void checkForcesAgainstReference(const TrajectoryFrame&              frame,
                                        const TrajectoryFrameMatchSettings& matchSettings,
                                        const TrajectoryTolerances&         tolerance,
                                        TestReferenceChecker*               checker)
{
    SCOPED_TRACE("Comparing forces");
    if (shouldDoComparison(frame.f(), matchSettings.forcesComparison))
    {
        checker->setDefaultTolerance(tolerance.forces);
        std::string id = frame.frameName() + " F";
        checker->checkSequence(std::begin(frame.f()), std::end(frame.f()), id.c_str());
    }
}

/*! \brief Compares the box from \c frame to reference data
 * according to the \c matchSettings and \c tolerance.
 */
static void checkBoxAgainstReference(const TrajectoryFrame&              frame,
                                     const TrajectoryFrameMatchSettings& matchSettings,
                                     const TrajectoryTolerances&         tolerance,
                                     TestReferenceChecker*               checker)
{
    SCOPED_TRACE("Comparing box");
    if (!matchSettings.mustCompareBox)
    {
        return;
    }
    if (!frame.hasBox())
    {
        ADD_FAILURE() << "Comparing the box was required, "
                         "but the test frame did not have one";
        return;
    }
    checker->setDefaultTolerance(tolerance.box);
    // Do the comparison
    for (int d = 0; d < DIM; ++d)
    {
        for (int dd = 0; dd < DIM; ++dd)
        {
            checker->checkReal(
                    frame.box()[d][dd],
                    (frame.frameName() + " Box[" + toString(d) + "][" + toString(dd) + "]").c_str());
        }
    }
}

void TrajectoryComparison::operator()(const TrajectoryFrame& frame, TestReferenceChecker* checker) const
{
    SCOPED_TRACE("Comparing trajectory frame " + frame.frameName() + " against reference data");
    checkBoxAgainstReference(frame, matchSettings_, tolerances_, checker);
    checkPositionsAgainstReference(frame, matchSettings_, tolerances_, checker);
    checkVelocitiesAgainstReference(frame, matchSettings_, tolerances_, checker);
    checkForcesAgainstReference(frame, matchSettings_, tolerances_, checker);
}

void checkTrajectoryAgainstReferenceData(const std::string&          trajectoryFilename,
                                         const TrajectoryComparison& trajectoryComparison,
                                         TestReferenceChecker*       checker)
{
    checkTrajectoryAgainstReferenceData(
            trajectoryFilename, trajectoryComparison, checker, MaxNumFrames::compareAllFrames());
}

void checkTrajectoryAgainstReferenceData(const std::string&          trajectoryFilename,
                                         const TrajectoryComparison& trajectoryComparison,
                                         TestReferenceChecker*       checker,
                                         MaxNumFrames                maxNumFrames)
{
    ;
    TrajectoryFrameReader reader(trajectoryFilename);
    unsigned int          counter = 0;
    for (counter = 0; counter < maxNumFrames && reader.readNextFrame(); counter++)
    {
        auto frame = reader.frame();
        trajectoryComparison(frame, checker);
    }

    if (counter == maxNumFrames && reader.readNextFrame())
    {
        // There would have been at least one more frame!
        checker->disableUnusedEntriesCheck();
    }
}

} // namespace test
} // namespace gmx
