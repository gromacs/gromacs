/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief Langevin tests runners.
 *
 * Declares test runner class for Langevin integrator algorithm. The test runners abstract
 * class is used to unify the interfaces for CPU and GPU implementations of the
 * Langevin algorithm. This allows to run the same test on the same data using
 * different implementations of the parent class, that inherit its interfaces.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Magnus Lundborg <magnus.lundborg@scilifelab.se>
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_TESTS_LANGEVINTESTRUNNERS_H
#define GMX_MDLIB_TESTS_LANGEVINTESTRUNNERS_H

#include "config.h"

#include <string>

#include "gromacs/math/vec.h"

#include "testutils/test_device.h"

#include "langevintestdata.h"

namespace gmx
{
namespace test
{

/* \brief Langevin integrator test runner interface.
 *
 * Wraps the actual implementation of Langevin algorithm into common interface.
 */
class ILangevinTestRunner
{
public:
    //! Virtual destructor
    virtual ~ILangevinTestRunner() {}
    /*! \brief The abstract function that runs the integrator for a given number of steps.
     *
     * Should be overriden.
     *
     * \param[in]     testData  Data needed for the integrator
     * \param[in]     numSteps  Total number of steps to run integration for.
     */
    virtual void integrate(LangevinTestData* testData, int numSteps) = 0;

    /*! \brief Get the human-friendly description of hardware used by the runner.
     *
     * \returns String with description of the hardware.
     */
    virtual std::string hardwareDescription() = 0;
};

// Runner for the CPU version of Langevin integrator.
class LangevinHostTestRunner : public ILangevinTestRunner
{
public:
    //! Constructor.
    LangevinHostTestRunner() {}
    /*! \brief Integrate on the CPU for a given number of steps.
     *
     * Will update the test data with the integration result.
     *
     * \param[in]     testData  Data needed for the integrator
     * \param[in]     numSteps  Total number of steps to run integration for.
     */
    void integrate(LangevinTestData* testData, int numSteps) override;
    /*! \brief Get the hardware description
     *
     * \returns "CPU" string.
     */
    std::string hardwareDescription() override { return "CPU"; }
};

} // namespace test
} // namespace gmx

#endif // GMX_MDLIB_TESTS_LANGEVINTESTRUNNERS_H
