/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
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
 * \brief Declares types and functions for comparing trajectories
 * produced by mdrun.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#ifndef GMX_PROGRAMS_MDRUN_TESTS_TRAJECTORYCOMPARISON_H
#define GMX_PROGRAMS_MDRUN_TESTS_TRAJECTORYCOMPARISON_H

#include "testutils/testasserts.h"

namespace gmx
{

class TrajectoryFrame;

namespace test
{

/*! \internal
 * \brief Helper struct for testing different trajectory components with different tolerances. */
struct TrajectoryTolerances
{
    /*!@{*/
    /*! \brief Tolerances for reproduction of different quantities. */
    FloatingPointTolerance box, coordinates, velocities, forces;
    /*!@}*/
};

//! Enumeration controlling how data within trajectory frames are compared
enum class ComparisonConditions : int
{
    CompareIfBothFound,
    NoComparison,
    MustCompare,
    CompareIfReferenceFound,
    CompareIfTestFound,
    Count
};

/*! \internal
 * \brief Helper struct to specify the expected behaviour of compareFrames().
 *
 * By default, nothing is required to be compared, but the comparer will
 * compare what it can with the frames it is given.
 *
 * Handling PBC refers to putting all the atoms in the simulation box,
 * which requires that both the PBC type and a simulation box are
 * available from the trajectory frame. */
struct TrajectoryFrameMatchSettings
{
    //! Whether boxes must be compared.
    bool mustCompareBox = false;
    //! Whether PBC will be handled if it can be handled.
    bool handlePbcIfPossible = true;
    //! Whether PBC handling must occur for a valid comparison.
    bool requirePbcHandling = false;
    //! Whether position coordinates must be compared.
    ComparisonConditions coordinatesComparison = ComparisonConditions::CompareIfBothFound;
    //! Whether velocities must be compared.
    ComparisonConditions velocitiesComparison = ComparisonConditions::CompareIfBothFound;
    //! Whether forces must be compared.
    ComparisonConditions forcesComparison = ComparisonConditions::CompareIfBothFound;
};

/*! \internal
 * \brief Function object to compare the fields of the two frames for
 * equality given the \c matchSettings_ and \c tolerances_.
 *
 * The two frames are required to have valid and matching values for
 * time and step. According to \c matchSettings_, box, position coordinates,
 * velocities and/or forces will be compared between frames, using the
 * \c tolerances_. Comparisons will only occur when both frames have
 * the requisite data, and will be expected to be equal within the
 * matching component of \c tolerances_. If a comparison fails, a
 * GoogleTest expectation failure will be given. If a comparison is
 * required by \c matchSettings_ but cannot be done because either (or
 * both) frames lack the requisite data, descriptive expectation
 * failures will be given. */
class TrajectoryComparison
{
public:
    //! Defaults for trajectory comparisons
    static const TrajectoryTolerances s_defaultTrajectoryTolerances;
    //! Constructor
    TrajectoryComparison(const TrajectoryFrameMatchSettings& matchSettings,
                         const TrajectoryTolerances&         tolerances);
    /*! \brief Compare reference with test given the \c
     * matchSettings_ within \c tolerances_ */
    void operator()(const TrajectoryFrame& reference, const TrajectoryFrame& test) const;
    //! Specifies expected behavior in comparisons
    TrajectoryFrameMatchSettings matchSettings_;
    //! Trajectory fields to match with given tolerances.
    TrajectoryTolerances tolerances_;
};

} // namespace test
} // namespace gmx

#endif
