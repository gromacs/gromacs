/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019,2020, by the GROMACS development team, led by
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
#ifndef GMX_TESTUTILS_TEST_HARDWARE_ENVIRONMENT_H
#define GMX_TESTUTILS_TEST_HARDWARE_ENVIRONMENT_H

/*! \internal \file
 * \brief
 * Describes test environment class which performs hardware enumeration for unit tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_testutils
 */

#include <map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/gmxassert.h"

#include "testutils/test_device.h"

struct gmx_hw_info_t;

namespace gmx
{
namespace test
{

/*! \internal \brief
 * This class performs one-time test initialization (enumerating the hardware)
 */
class TestHardwareEnvironment : public ::testing::Environment
{
private:
    //! General hardware info
    gmx_hw_info_t* hardwareInfo_;
    //! Storage of hardware contexts
    std::vector<std::unique_ptr<TestDevice>> testDeviceList_;

public:
    //! This is called by GTest framework once to query the hardware
    void SetUp() override;
    //! This is called by GTest framework once release the hardware
    void TearDown() override;
    //! Get available hardware contexts.
    const std::vector<std::unique_ptr<TestDevice>>& getTestDeviceList() const
    {
        return testDeviceList_;
    }
    bool hasCompatibleDevices() const { return !testDeviceList_.empty(); }
    //! Get available hardware information.
    const gmx_hw_info_t* hwinfo() const { return hardwareInfo_; }
};

//! Get the test environment
const TestHardwareEnvironment* getTestHardwareEnvironment();

/*! \brief This constructs the test environment during setup of the
 * unit test so that they can use the hardware context. */
void callAddGlobalTestEnvironment();

} // namespace test
} // namespace gmx
#endif // GMX_TESTUTILS_TEST_HARDWARE_ENVIRONMENT_H
