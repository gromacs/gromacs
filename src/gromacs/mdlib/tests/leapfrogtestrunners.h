/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * \brief Leap-Frog tests runners.
 *
 * Declares test runner class for Leap-Frog algorithm. The test runners abstract
 * class is used to unify the interfaces for CPU and GPU implementations of the
 * Leap-Frog algorithm. This allows to run the same test on the same data using
 * different implementations of the parent class, that inherit its interfaces.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_TESTS_LEAPFROGTESTRUNNERS_H
#define GMX_MDLIB_TESTS_LEAPFROGTESTRUNNERS_H

#include "config.h"

#include <string>

#include "gromacs/math/vec.h"

#include "testutils/test_device.h"

#include "leapfrogtestdata.h"

/*
 * LeapFrog is available with CUDA and SYCL.
 */
#define GPU_LEAPFROG_SUPPORTED (GMX_GPU_CUDA || GMX_GPU_SYCL)

namespace gmx
{
namespace test
{

/* \brief LeapFrog integrator test runner interface.
 *
 * Wraps the actual implementation of LeapFrog algorithm into common interface.
 */
class ILeapFrogTestRunner
{
public:
    //! Virtual destructor
    virtual ~ILeapFrogTestRunner() {}
    /*! \brief The abstract function that runs the integrator for a given number of steps.
     *
     * Should be overriden.
     *
     * \param[in]     testData  Data needed for the integrator
     * \param[in]     numSteps  Total number of steps to run integration for.
     */
    virtual void integrate(LeapFrogTestData* testData, int numSteps) = 0;

    /*! \brief Get the human-friendly description of hardware used by the runner.
     *
     * \returns String with description of the hardware.
     */
    virtual std::string hardwareDescription() = 0;
};

// Runner for the CPU version of Leap-Frog.
class LeapFrogHostTestRunner : public ILeapFrogTestRunner
{
public:
    //! Constructor.
    LeapFrogHostTestRunner() {}
    /*! \brief Integrate on the CPU for a given number of steps.
     *
     * Will update the test data with the integration result.
     *
     * \param[in]     testData  Data needed for the integrator
     * \param[in]     numSteps  Total number of steps to run integration for.
     */
    void integrate(LeapFrogTestData* testData, int numSteps) override;
    /*! \brief Get the hardware description
     *
     * \returns "CPU" string.
     */
    std::string hardwareDescription() override { return "CPU"; }
};

// Runner for the CPU version of Leap-Frog.
class LeapFrogDeviceTestRunner : public ILeapFrogTestRunner
{
public:
    /*! \brief Constructor. Keeps a copy of the hardware context.
     *
     * \param[in] testDevice The device hardware context to be used by the runner.
     */
    LeapFrogDeviceTestRunner(const TestDevice& testDevice) : testDevice_(testDevice) {}
    /*! \brief Integrate on the GPU for a given number of steps.
     *
     * Copies data from CPU to GPU, integrates the equation of motion
     * for requested number of steps using Leap-Frog algorithm, copies
     * the result back.
     *
     * \param[in]     testData  Data needed for the integrator
     * \param[in]     numSteps  Total number of steps to run integration for.
     */
    void integrate(LeapFrogTestData* testData, int numSteps) override;
    /*! \brief Get the hardware description
     *
     * \returns A string with GPU description.
     */
    std::string hardwareDescription() override { return testDevice_.description(); }

private:
    //! Test device to be used in the runner.
    const TestDevice& testDevice_;
};

} // namespace test
} // namespace gmx

#endif // GMX_MDLIB_TESTS_LEAPFROGTESTRUNNERS_H
