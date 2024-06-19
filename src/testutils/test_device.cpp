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
 * \brief
 * Implements test environment class which performs hardware enumeration for unit tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/test_device.h"

#include <memory>
#include <string>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/device_information.h"

namespace gmx
{
namespace test
{

class TestDevice::Impl
{
public:
    Impl(const char* description, const DeviceInformation& deviceInfo);
    ~Impl();
    //! Returns a human-readable context description line
    std::string description() const { return description_; }
    //! Returns the device info pointer
    const DeviceInformation& deviceInfo() const { return deviceContext_.deviceInfo(); }
    //! Get the device context
    const DeviceContext& deviceContext() const { return deviceContext_; }
    //! Get the device stream
    const DeviceStream& deviceStream() const { return deviceStream_; }

private:
    //! Readable description
    std::string description_;
    //! Device context
    DeviceContext deviceContext_;
    //! Device stream
    DeviceStream deviceStream_;
};

TestDevice::Impl::Impl(const char* description, const DeviceInformation& deviceInfo) :
    description_(description),
    deviceContext_(deviceInfo),
    deviceStream_(deviceContext_, DeviceStreamPriority::Normal, false)
{
}

TestDevice::Impl::~Impl() = default;

TestDevice::TestDevice(const char* description, const DeviceInformation& deviceInfo) :
    impl_(new Impl(description, deviceInfo))
{
}

TestDevice::~TestDevice() = default;

std::string TestDevice::description() const
{
    return impl_->description();
}

int TestDevice::id() const
{
    return deviceInfo().id;
}

const DeviceInformation& TestDevice::deviceInfo() const
{
    return impl_->deviceInfo();
}

const DeviceContext& TestDevice::deviceContext() const
{
    return impl_->deviceContext();
}

const DeviceStream& TestDevice::deviceStream() const
{
    return impl_->deviceStream();
}

} // namespace test
} // namespace gmx
