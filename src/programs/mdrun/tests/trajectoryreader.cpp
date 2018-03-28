/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018, by the GROMACS development team, led by
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
 * \brief Implementions of related classes for tests that want to
 * inspect trajectories produced by mdrun.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "trajectoryreader.h"

#include <memory>
#include <string>

#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include <gmock/gmock.h>

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

/*! \brief Implementation class for RealEq matcher
 *
 * See RealEq().
 */
template <typename FloatType>
class FloatTypeMatcher : public testing::MatcherInterface < std::tuple < FloatType, FloatType>>
{
    public:
        //! Constructor
        FloatTypeMatcher(FloatingPointTolerance tolerance)
            : tolerance_(tolerance) {}
        //! Compare the two elements of \c arg, return whether they are equal, and comment on \c listener when they are not.
        virtual bool MatchAndExplain(std::tuple<FloatType, FloatType> arg,
                                     testing::MatchResultListener* listener) const
        {
            const FloatType        &value1 = std::get<0>(arg);
            const FloatType        &value2 = std::get<1>(arg);
            FloatingPointDifference diff(value1, value2);
            if (tolerance_.isWithin(diff))
            {
                return true;
            }
            *listener->stream()
            << "  Actual value: " << value2 << std::endl
            << "Expected value: " << value1 << std::endl
            << "    Difference: " << diff.toString() << std::endl
            << "     Tolerance: " << tolerance_.toString(diff);
            return false;
        }
        //! Describe to a human what matching means.
        virtual void DescribeTo(::std::ostream* os) const
        {
            *os << "matches within tolerance";
        }
        //! Describe to a human what failing to match means.
        virtual void DescribeNegationTo(::std::ostream* os) const
        {
            *os << "does not match within tolerance";
        }
    private:
        //! Tolerance used in matching
        FloatingPointTolerance tolerance_;
};

/*! \brief Make matcher for reals for use with GoogleMock that compare
 * equal when \c tolerance is satisifed.
 *
 * Used like
 *
 *   EXPECT_THAT(testReals, Pointwise(RealEq(tolerance), referenceReals));
 */
template <typename FloatType>
inline testing::Matcher < std::tuple < FloatType, FloatType>>
RealEq(FloatingPointTolerance tolerance)
{
    return testing::MakeMatcher(new FloatTypeMatcher<FloatType>(tolerance));
}

/*! \brief Implementation class for RvecEq matcher
 *
 * See RvecEq().
 */
template <typename FloatType>
class RVecMatcher :
    public testing::MatcherInterface < std::tuple < BasicVector<FloatType>, BasicVector<FloatType>>>
{
    public:
        //! Convenience type
        using VectorType = BasicVector<FloatType>;
        //! Constructor
        RVecMatcher(FloatingPointTolerance tolerance)
            : tolerance_(tolerance) {}
        //! Compare the two elements of \c arg, return whether they are equal, and comment on \c listener when they are not.
        virtual bool MatchAndExplain(std::tuple<VectorType, VectorType> arg,
                                     testing::MatchResultListener* listener) const
        {
            const VectorType           &lhs = std::get<0>(arg);
            const VectorType           &rhs = std::get<1>(arg);
            FloatTypeMatcher<FloatType> floatTypeMatcher(tolerance_);
            bool matches = true;
            for (int d = 0; d < DIM; ++d)
            {
                auto floatTuple = std::make_tuple<FloatType, FloatType>(lhs[d], rhs[d]);
                matches = matches && floatTypeMatcher.MatchAndExplain(floatTuple, listener);
            }
            return matches;
        }
        //! Describe to a human what matching means.
        virtual void DescribeTo(::std::ostream* os) const
        {
            *os << "matches all elements within tolerance";
        }
        //! Describe to a human what failing to match means.
        virtual void DescribeNegationTo(::std::ostream* os) const
        {
            *os << "does not match all elements within tolerance";
        }
    private:
        //! Tolerance used in matching
        FloatingPointTolerance tolerance_;
};

/*! \brief Make matcher for RVecs for use with GoogleMock that compare
 * equal when \c tolerance is satisifed.
 *
 * Used like
 *
 *   EXPECT_THAT(testRVecs, Pointwise(RVecEq(tolerance), referenceRVecs));
 */
inline testing::Matcher < std::tuple < RVec, RVec>>
RVecEq(FloatingPointTolerance tolerance)
{
    return testing::MakeMatcher(new RVecMatcher<real>(tolerance));
}

//! Helper function to obtain resources
static t_trxframe *make_trxframe()
{
    t_trxframe *frame;

    snew(frame, 1);
    clear_trxframe(frame, true);

    return frame;
}

//! Helper function to clean up resources
void done_trxframe(t_trxframe *fr)
{
    // Free the contents, then the pointer itself
    sfree(fr->x);
    sfree(fr->v);
    sfree(fr->f);
    sfree(fr->index);
    sfree(fr);
}

// === TrajectoryFrameReader ===

TrajectoryFrameReader::TrajectoryFrameReader(const std::string &filename)
    : filename_(filename),
      trajectoryFileGuard_(nullptr),
      trxframeGuard_(make_trxframe()),
      haveReadFirstFrame_(false),
      haveProbedForNextFrame_(false),
      nextFrameExists_(false)
{
    gmx_output_env_t *oenv;
    output_env_init_default(&oenv);
    oenvGuard_.reset(oenv);
}

bool
TrajectoryFrameReader::readNextFrame()
{
    if (haveProbedForNextFrame_)
    {
        if (nextFrameExists_)
        {
            GMX_THROW(APIError("This frame has already been probed for, it should be used before probing again."));
        }
        else
        {
            GMX_THROW(APIError("This frame has already been probed for, it doesn't exist, so there should not be subsequent attempts to probe for it."));
        }
    }
    haveProbedForNextFrame_ = true;
    // If there's a next frame, read it into trxframe_, and report the result.
    if (!haveReadFirstFrame_)
    {
        t_trxstatus *trajectoryFile;
        int          flags = TRX_READ_X | TRX_READ_V | TRX_READ_F;
        nextFrameExists_ = read_first_frame(oenvGuard_.get(),
                                            &trajectoryFile,
                                            filename_.c_str(),
                                            trxframeGuard_.get(),
                                            flags);
        if (!trajectoryFile)
        {
            GMX_THROW(FileIOError("Could not open trajectory file " + filename_ + " for reading"));
        }
        trajectoryFileGuard_.reset(trajectoryFile);
        haveReadFirstFrame_ = true;
    }
    else
    {
        nextFrameExists_ = read_next_frame(oenvGuard_.get(),
                                           trajectoryFileGuard_.get(),
                                           trxframeGuard_.get());
    }
    return nextFrameExists_;
}

TrajectoryFrame
TrajectoryFrameReader::frame()
{
    if (!haveProbedForNextFrame_)
    {
        readNextFrame();
    }
    if (!nextFrameExists_)
    {
        GMX_THROW(APIError("There is no next frame, so there should have been no attempt to get it. Perhaps the return value of readNextFrame() was misused."));
    }

    // Prepare for reading future frames
    haveProbedForNextFrame_ = false;
    nextFrameExists_        = false;

    // The probe filled trxframeGuard_ with new data, so return it
    return TrajectoryFrame(*trxframeGuard_.get());
}

// === Free functions ===

using ::testing::Pointwise;

/*! \brief Compares the box from \c reference and \c test
 * according to the \c matchSettings and \c tolerance.
 *
 * \todo This could be streamlined when we have a proper 3D matrix
 * class and view. */
static void compareBox(const TrajectoryFrame              &reference,
                       const TrajectoryFrame              &test,
                       const TrajectoryFrameMatchSettings &matchSettings,
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

/*! \brief Help put all atom positions in \c frame into its box.
 *
 * This can perhaps go away when frame->x is a container. */
static std::vector<RVec>
putAtomsInBox(const TrajectoryFrame &frame)
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

/*! \brief Compares the positions from \c reference and \c test
 * according to the \c matchSettings and \c tolerance. */
static void comparePositions(const TrajectoryFrame              &reference,
                             const TrajectoryFrame              &test,
                             const TrajectoryFrameMatchSettings &matchSettings,
                             const FloatingPointTolerance        tolerance)
{
    bool canHandlePbc = true;
    if (!reference.hasBox())
    {
        if (matchSettings.mustComparePositions)
        {
            ADD_FAILURE() << "Comparing positions required PBC handling, "
            "but the reference frame did not have a box";
        }
        canHandlePbc = false;
    }
    if (!test.hasBox())
    {
        if (matchSettings.mustComparePositions)
        {
            ADD_FAILURE() << "Comparing positions required PBC handling, "
            "but the test frame did not have a box";
        }
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
static void compareVelocities(const TrajectoryFrame              &reference,
                              const TrajectoryFrame              &test,
                              const TrajectoryFrameMatchSettings &matchSettings,
                              const FloatingPointTolerance        tolerance)
{
    if (!matchSettings.mustCompareVelocities)
    {
        return;
    }
    EXPECT_THAT(test.v(), Pointwise(RVecEq(tolerance), reference.v()));
}

/*! \brief Compares the forces from \c reference and \c test
 * according to the \c matchSettings and \c tolerance. */
static void compareForces(const TrajectoryFrame              &reference,
                          const TrajectoryFrame              &test,
                          const TrajectoryFrameMatchSettings &matchSettings,
                          const FloatingPointTolerance        tolerance)
{
    if (!matchSettings.mustCompareForces)
    {
        return;
    }
    EXPECT_THAT(test.f(), Pointwise(RVecEq(tolerance), reference.f()));
}


void compareTrajectoryFrames(const TrajectoryFrame              &reference,
                             const TrajectoryFrame              &test,
                             const TrajectoryFrameMatchSettings &matchSettings,
                             const TrajectoryTolerances         &tolerances)
{
    SCOPED_TRACE("Comparing reference frame " + reference.frameName() + " and test frame " + test.frameName());
    EXPECT_EQ(reference.step(), test.step());
    EXPECT_EQ(reference.time(), test.time());
    compareBox(reference, test, matchSettings, tolerances.box);
    comparePositions(reference, test, matchSettings, tolerances.positions);
    compareVelocities(reference, test, matchSettings, tolerances.velocities);
    compareForces(reference, test, matchSettings, tolerances.forces);
}

} // namespace
} // namespace
