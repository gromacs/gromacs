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
 * \brief SHAKE and LINCS tests runners.
 *
 * Declares test runner class for constraints. The test runner abstract class is used
 * to unify the interfaces for different constraints methods, running on different
 * hardware.  This allows to run the same test on the same data using different
 * implementations of the parent class, that inherit its interfaces.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_TESTS_CONSTRTESTRUNNERS_H
#define GMX_MDLIB_TESTS_CONSTRTESTRUNNERS_H

#include <string>

#include <gtest/gtest.h>

#include "gromacs/pbcutil/pbc.h"

#include "testutils/test_device.h"

#include "constrtestdata.h"

/*
 * GPU version of constraints is only available with CUDA and SYCL.
 */
#define GPU_CONSTRAINTS_SUPPORTED (GMX_GPU_CUDA || GMX_GPU_SYCL)

struct t_pbc;

namespace gmx
{
namespace test
{

/* \brief Constraints test runner interface.
 *
 * Wraps the actual implementation of constraints algorithm into common interface.
 */
class IConstraintsTestRunner
{
public:
    //! Virtual destructor.
    virtual ~IConstraintsTestRunner() {}
    /*! \brief Abstract constraining function. Should be overriden.
     *
     * \param[in] testData             Test data structure.
     * \param[in] pbc                  Periodic boundary data.
     */
    virtual void applyConstraints(ConstraintsTestData* testData, t_pbc pbc) = 0;

    /*! \brief Get the name of the implementation.
     *
     * \return "<algorithm> on <device>", depending on the actual implementation used. E.g., "LINCS on #0: NVIDIA GeForce GTX 1660 SUPER".
     */
    virtual std::string name() = 0;
};

// Runner for the CPU implementation of SHAKE constraints algorithm.
class ShakeConstraintsRunner : public IConstraintsTestRunner
{
public:
    //! Default constructor.
    ShakeConstraintsRunner() {}
    /*! \brief Apply SHAKE constraints to the test data.
     *
     * \param[in] testData             Test data structure.
     * \param[in] pbc                  Periodic boundary data.
     */
    void applyConstraints(ConstraintsTestData* testData, t_pbc pbc) override;
    /*! \brief Get the name of the implementation.
     *
     * \return "SHAKE" string;
     */
    std::string name() override { return "SHAKE on CPU"; }
};

// Runner for the CPU implementation of LINCS constraints algorithm.
class LincsConstraintsRunner : public IConstraintsTestRunner
{
public:
    //! Default constructor.
    LincsConstraintsRunner() {}
    /*! \brief Apply LINCS constraints to the test data on the CPU.
     *
     * \param[in] testData             Test data structure.
     * \param[in] pbc                  Periodic boundary data.
     */
    void applyConstraints(ConstraintsTestData* testData, t_pbc pbc) override;
    /*! \brief Get the name of the implementation.
     *
     * \return "LINCS" string;
     */
    std::string name() override { return "LINCS on CPU"; }
};

// Runner for the GPU implementation of LINCS constraints algorithm.
class LincsDeviceConstraintsRunner : public IConstraintsTestRunner
{
public:
    /*! \brief Constructor. Keeps a copy of the hardware context.
     *
     * \param[in] testDevice The device hardware context to be used by the runner.
     */
    LincsDeviceConstraintsRunner(const TestDevice& testDevice) : testDevice_(testDevice) {}
    /*! \brief Apply LINCS constraints to the test data on the GPU.
     *
     * \param[in] testData             Test data structure.
     * \param[in] pbc                  Periodic boundary data.
     */
    void applyConstraints(ConstraintsTestData* testData, t_pbc pbc) override;
    /*! \brief Get the name of the implementation.
     *
     * \return "LINCS_GPU" string;
     */
    std::string name() override { return "LINCS on " + testDevice_.description(); }

private:
    //! Test device to be used in the runner.
    const TestDevice& testDevice_;
};

} // namespace test
} // namespace gmx

#endif // GMX_MDLIB_TESTS_CONSTRTESTRUNNERS_H
