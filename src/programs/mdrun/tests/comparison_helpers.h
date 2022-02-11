/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * \brief Declares types and functions common to comparing either
 * energies or trajectories produced by mdrun.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#ifndef GMX_PROGRAMS_MDRUN_TESTS_COMPARISON_HELPERS_H
#define GMX_PROGRAMS_MDRUN_TESTS_COMPARISON_HELPERS_H

#include <limits>

namespace gmx
{

namespace test
{

/*! \internal
 * \brief Named struct indicating the max number of frames to be compared */
struct MaxNumFrames
{
    //! Explicit constructor
    explicit MaxNumFrames(unsigned int maxFrame) : maxFrame_(maxFrame) {}

    //! Implicit conversion to int - struct can be used like underlying type
    operator unsigned int() const { return maxFrame_; }

    //! Return a MaxNumFrames that will try to compare all frames
    [[nodiscard]] static MaxNumFrames compareAllFrames()
    {
        return MaxNumFrames(std::numeric_limits<decltype(maxFrame_)>::max());
    }

private:
    //! Internal value
    const unsigned int maxFrame_;
};

} // namespace test
} // namespace gmx

#endif
