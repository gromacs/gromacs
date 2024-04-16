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
#ifndef GMX_TESTUTILS_TEST_DEVICE_H
#define GMX_TESTUTILS_TEST_DEVICE_H

/*! \internal \file
 * \brief
 * Describes test environment class which performs GPU device enumeration for unit tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_testutils
 */

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/utility/gmxassert.h"

class DeviceContext;
struct DeviceInformation;
class DeviceStream;

namespace gmx
{
namespace test
{

/*! \internal \brief
 * A structure to describe a hardware context that persists over the lifetime
 * of the test binary.
 */
class TestDevice
{
public:
    //! Returns a human-readable context description line
    std::string description() const;
    //! Returns a numerical ID for the device
    int id() const;
    //! Returns the device info pointer
    const DeviceInformation& deviceInfo() const;
    //! Get the device context, particularly useful to activate the device
    const DeviceContext& deviceContext() const;
    //! Get the device stream
    const DeviceStream& deviceStream() const;
    //! Creates the device context and stream for tests on the GPU
    TestDevice(const char* description, const DeviceInformation& deviceInfo);
    //! Destructor
    ~TestDevice();

private:
    //! Implementation type.
    class Impl;
    //! Implementation object.
    std::unique_ptr<Impl> impl_;
};

} // namespace test
} // namespace gmx

#endif // GMX_TESTUTILS_TEST_DEVICE_H
