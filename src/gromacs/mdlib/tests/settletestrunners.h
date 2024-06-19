/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief SETTLE tests runners.
 *
 * Declares test runner class for SETTLE algorithm. The test runners abstract
 * class is used to unify the interfaces for CPU and GPU implementations of the
 * SETTLE algorithm. This allows to run the same test on the same data using
 * different implementations of the parent class, that inherit its interfaces.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_TESTS_SETTLETESTRUNNERS_H
#define GMX_MDLIB_TESTS_SETTLETESTRUNNERS_H

#include <string>

#include <gtest/gtest.h>

#include "gromacs/pbcutil/pbc.h"

#include "testutils/test_device.h"

#include "settletestdata.h"

/*
 * GPU version of SETTLE is only available with CUDA.
 */
#define GPU_SETTLE_SUPPORTED (GMX_GPU_CUDA || GMX_GPU_SYCL)

struct t_pbc;

namespace gmx
{
namespace test
{

/* \brief SETTLE test runner interface.
 *
 * Wraps the actual implementation of SETTLE into common interface.
 */
class ISettleTestRunner
{
public:
    //! Virtual destructor.
    virtual ~ISettleTestRunner() {}

    /*! \brief Apply SETTLE using CPU version of the algorithm
     *
     * Initializes SETTLE object, applies algorithm, destroys the object. The coordinates, velocities
     * and virial are updated in the testData object.
     *
     * \param[in,out] testData          An object, containing all the data structures needed by SETTLE.
     * \param[in]     pbc               Periodic boundary setup.
     * \param[in]     updateVelocities  If the velocities should be updated.
     * \param[in]     calcVirial        If the virial should be computed.
     * \param[in]     testDescription   Brief description that will be printed in case of test failure.
     */
    virtual void applySettle(SettleTestData*    testData,
                             t_pbc              pbc,
                             bool               updateVelocities,
                             bool               calcVirial,
                             const std::string& testDescription) = 0;
    /*! \brief Get the hardware description
     *
     * \returns A string, describing hardware used by the runner.
     */
    virtual std::string hardwareDescription() = 0;
};

// Runner for the CPU implementation of SETTLE.
class SettleHostTestRunner : public ISettleTestRunner
{
public:
    //! Default constructor.
    SettleHostTestRunner() {}
    /*! \brief Apply SETTLE using CPU version of the algorithm
     *
     * Initializes SETTLE object, applies algorithm, destroys the object. The coordinates, velocities
     * and virial are updated in the testData object.
     *
     * \param[in,out] testData          An object, containing all the data structures needed by SETTLE.
     * \param[in]     pbc               Periodic boundary setup.
     * \param[in]     updateVelocities  If the velocities should be updated.
     * \param[in]     calcVirial        If the virial should be computed.
     * \param[in]     testDescription   Brief description that will be printed in case of test failure.
     */
    void applySettle(SettleTestData*    testData,
                     t_pbc              pbc,
                     bool               updateVelocities,
                     bool               calcVirial,
                     const std::string& testDescription) override;
    /*! \brief Get the hardware description
     *
     * \returns "CPU" string.
     */
    std::string hardwareDescription() override { return "CPU"; }
};

// Runner for the GPU implementation of SETTLE.
class SettleDeviceTestRunner : public ISettleTestRunner
{
public:
    /*! \brief Constructor. Keeps a copy of the hardware context.
     *
     * \param[in] testDevice The device hardware context to be used by the runner.
     */
    SettleDeviceTestRunner(const TestDevice& testDevice) : testDevice_(testDevice) {}
    /*! \brief Apply SETTLE using GPU version of the algorithm
     *
     * Initializes SETTLE object, copied data to the GPU, applies algorithm, copies the data back,
     * destroys the object. The coordinates, velocities and virial are updated in the testData object.
     *
     * \param[in,out] testData          An object, containing all the data structures needed by SETTLE.
     * \param[in]     pbc               Periodic boundary setup.
     * \param[in]     updateVelocities  If the velocities should be updated.
     * \param[in]     calcVirial        If the virial should be computed.
     * \param[in]     testDescription   Brief description that will be printed in case of test failure.
     */
    void applySettle(SettleTestData*    testData,
                     t_pbc              pbc,
                     bool               updateVelocities,
                     bool               calcVirial,
                     const std::string& testDescription) override;
    /*! \brief Get the hardware description
     *
     * \returns A string with GPU description.
     */
    std::string hardwareDescription() override { return testDevice_.description(); }

private:
    //! Test test device to be used in the runner.
    const TestDevice& testDevice_;
};

} // namespace test
} // namespace gmx

#endif // GMX_MDLIB_TESTS_SETTLETESTRUNNERS_H
