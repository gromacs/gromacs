/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements test environment class which performs hardware enumeration for unit tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include "testhardwarecontext.h"

#include <memory>

#include "gromacs/ewald/pme.h"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/loggerbuilder.h"
#include "gromacs/utility/physicalnodecommunicator.h"

namespace gmx
{
namespace test
{

TestHardwareContext::TestHardwareContext(CodePath codePath, const char* description) :
    codePath_(codePath),
    description_(description)
{
    GMX_RELEASE_ASSERT(codePath == CodePath::CPU,
                       "A GPU code path should provide DeviceInformation to the "
                       "TestHerdwareContext constructor.");
    deviceContext_ = nullptr;
    deviceStream_  = nullptr;
}

TestHardwareContext::TestHardwareContext(CodePath                 codePath,
                                         const char*              description,
                                         const DeviceInformation& deviceInfo) :
    codePath_(codePath),
    description_(description)
{
    GMX_RELEASE_ASSERT(codePath == CodePath::GPU,
                       "TestHardwareContext tries to construct DeviceContext and PmeGpuProgram "
                       "in CPU build.");
    deviceContext_ = new DeviceContext(deviceInfo);
    deviceStream_  = new DeviceStream(*deviceContext_, DeviceStreamPriority::Normal, false);
    program_       = buildPmeGpuProgram(*deviceContext_);
}

TestHardwareContext::~TestHardwareContext()
{
    delete (deviceStream_);
    delete (deviceContext_);
}

const DeviceInformation* TestHardwareContext::deviceInfo() const
{
    return &deviceContext_->deviceInfo();
}

const DeviceContext* TestHardwareContext::deviceContext() const
{
    return deviceContext_;
}
//! Get the device stream
const DeviceStream* TestHardwareContext::deviceStream() const
{
    return deviceStream_;
}

const char* codePathToString(CodePath codePath)
{
    switch (codePath)
    {
        case CodePath::CPU: return "CPU";
        case CodePath::GPU: return "GPU";
        default: GMX_THROW(NotImplementedError("This CodePath should support codePathToString"));
    }
}

} // namespace test
} // namespace gmx
