/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \libinternal
 *
 * \brief Functionality for testing whether calls to mdrun produce the
 * same energy quantities when they should do so.
 */
#ifndef GMX_MDRUN_TESTS_MDRUNENERGYCOMPARISONFIXTURE_H
#define GMX_MDRUN_TESTS_MDRUNENERGYCOMPARISONFIXTURE_H

#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/integrationtests.h"

#include "moduletest.h"

namespace gmx
{

namespace test
{

class EnergyFrameInfo;
class FloatingPointTolerance;

/*! \libinternal \brief Declares test fixture class for integration
 * tests of mdrun functionality that will compare energies produced by
 * multiple calls to mdrun.
 *
 * A database of several kinds of simulation useful for such
 * comparisons is available.
 *
 * Any method in this class may throw std::bad_alloc if out of memory.
 *
 * \ingroup module_mdrun_integration_tests
 */
class MdrunEnergyComparisonFixture : public MdrunTestFixture
{
    public:
        //! Destructor
        virtual ~MdrunEnergyComparisonFixture();
        //! Helper typedef
        typedef std::pair<EnergyFrameInfo, EnergyFrameInfo> EnergyFrameInfoPair;
        /*! \brief Call grompp to do a simulation
         *
         * A database of several kinds of simulation useful for such
         * comparisons is available.
         *     - argon
         *     - spc216
         *
         * This method uses files in the database that match \c
         * simulationName, along with a matching .mdp string,
         * interpolating the other parameters into the fields of the
         * .mdp file.
         */
        void prepareSimulation(const char *simulationName,
                               const char *integrator,
                               const char *tcoupl,
                               const char *pcoupl);
        //! Compare the named fields in the two frames for equality within the tolerance
        void compareFrames(const EnergyFrameInfoPair    &frames,
                           std::vector<std::string>      namesOfRequiredFields,
                           const FloatingPointTolerance &tolerance);
};

} // namespace test
} // namespace gmx

#endif
