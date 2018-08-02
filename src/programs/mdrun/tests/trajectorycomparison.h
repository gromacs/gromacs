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
    FloatingPointTolerance box, positions, velocities, forces;
    /*!@}*/
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
    bool mustCompareBox;
    //! Whether positions must be compared.
    bool mustComparePositions;
    //! Whether PBC will be handled if it can be handled.
    bool handlePbcIfPossible;
    //! Whether PBC handling must occur for a valid comparison.
    bool requirePbcHandling;
    //! Whether velocities must be compared.
    bool mustCompareVelocities;
    //! Whether forces must be compared.
    bool mustCompareForces;
};

/*! \brief Compare the fields of the two frames for equality given
 * the \c matchSettings and \c tolerances.
 *
 * The two frames are required to have valid and matching values for
 * time and step. According to \c matchSettings, box, positions,
 * velocities and/or forces will be compared between frames, using the
 * \c tolerances. Comparisons will only occur when both frames have
 * the requisite data, and will be expected to be equal within the
 * matching component of \c tolerances. If a comparison fails, a
 * GoogleTest expectation failure will be given. If a comparison is
 * required by \c matchSettings but cannot be done because either (or
 * both) frames lack the requisite data, descriptive expectation
 * failures will be given. */
void compareTrajectoryFrames(const TrajectoryFrame              &reference,
                             const TrajectoryFrame              &test,
                             const TrajectoryFrameMatchSettings &matchSettings,
                             const TrajectoryTolerances         &tolerances);

}  // namespace test
}  // namespace gmx

#endif
