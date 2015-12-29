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
 * same energy and force quantities when they should do so.
 */
#ifndef GMX_MDRUN_TESTS_MDRUNCOMPARISONFIXTURE_H
#define GMX_MDRUN_TESTS_MDRUNCOMPARISONFIXTURE_H

#include <map>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"

#include "moduletest.h"

namespace gmx
{

namespace test
{

class EnergyFrame;
class TrajectoryFrame;
class FloatingPointTolerance;

/*! \libinternal
 * \brief Declares abstract base text fixture class for integration
 * tests of mdrun functionality that will compare multiple calls to
 * mdrun.
 *
 * An internal database of several kinds of simulation useful for such
 * comparisons is available.
 *
 * Any method in this class may throw std::bad_alloc if out of memory.
 *
 * \ingroup module_mdrun_integration_tests
 */
class MdrunComparisonFixture : public MdrunTestFixture
{
    public:
        //! Destructor
        virtual ~MdrunComparisonFixture();
        //! Helper typedef
        typedef std::map<std::string, std::string> MdpFieldValues;
        /*! \brief Prepare .mdp values to do a simulation
         *
         * A database of several kinds of simulation useful for different
         * kinds of tests is available.
         *     - argon12
         *     - argon5832
         *     - spc5
         *     - spc216
         *     - alanine_vsite_vacuo
         *     - alanine_vsite_solvated
         *     - nonanol
         *
         * Some of these systems are pretty minimal, because having
         * few atoms means few interactions, highly reproducible
         * forces, and allows tests to focus on the correctness of the
         * implementation of high-level mdrun features. The boxes are
         * of a reasonable size so that domain decomposition is
         * possible. The pressure-coupling parameters are isotropic,
         * and set up so that there will not be dramatic collapse of
         * volume over the handful of MD steps that will be run. A
         * single temperature-coupling group is used.
         *
         * This is separate from prepareMdpFile, so that derived
         * classes can react to the .mdp settings, e.g. by stopping a
         * run after half the steps.
         *
         * \throws  std::bad_alloc     if out of memory
         *          std::out_of_range  if \c simulationName is not in the database */
        MdpFieldValues prepareMdpFieldValues(const char *simulationName);
        /*! \brief Set up an .mdp file that permits a highly reproducible
         * simulation.
         *
         * \throws  std::bad_alloc     if out of memory */
        void prepareMdpFile(const MdpFieldValues &mdpFieldValues,
                            const char           *integrator,
                            const char           *tcoupl,
                            const char           *pcoupl);
        /*! \brief Run mdrun two ways in a test. Subclasses must override this method.
         *
         * It is expected that this method calls
         * prepareMdpFieldValues() and prepareMdpFile() to help set up
         * a call to grompp with gromppCallerRef. Then mdrun will be
         * called and perhaps energies and forces compared. */
        virtual void runTest(const CommandLine     &gromppCallerRef,
                             const char            *simulationName,
                             const char            *integrator,
                             const char            *tcoupl,
                             const char            *pcoupl,
                             FloatingPointTolerance tolerance) = 0;
        //! Convenience overload of runTest() for cases that don't need to customize the command line for grompp
        virtual void runTest(const char            *simulationName,
                             const char            *integrator,
                             const char            *tcoupl,
                             const char            *pcoupl,
                             FloatingPointTolerance tolerance);
};

} // namespace test
} // namespace gmx

#endif
