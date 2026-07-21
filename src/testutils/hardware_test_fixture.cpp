/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief Implements hardware-test fixture and utilities
 *
 * \author Mark Abraham
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/hardware_test_fixture.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/test_device.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

std::string TestHardwareContext::description() const
{
    return testDevice_ ? "GPU (" + testDevice_->description() + ")" : "CPU";
}

std::string TestHardwareContext::testName() const
{
    return testDevice_ ? "GPU" + std::to_string(testDevice_->id()) : "CPU";
}

std::optional<int> TestHardwareContext::gpuId() const
{
    return testDevice_ ? std::optional<int>(testDevice_->id()) : std::nullopt;
}

const DeviceContext* TestHardwareContext::deviceContext() const
{
    return testDevice_ ? &testDevice_->deviceContext() : nullptr;
}

const DeviceStream* TestHardwareContext::deviceStream() const
{
    return testDevice_ ? &testDevice_->deviceStream() : nullptr;
}

void TestHardwareContext::activate() const
{
    if (testDevice_)
    {
        testDevice_->deviceContext().activate();
    }
}

ArrayRef<const TestHardwareContext> getTestHardwareContexts()
{
    // Static storage with lazy initialization (avoid static initialization fiasco)
    static std::vector<TestHardwareContext> s_testHardwareContexts;

    if (s_testHardwareContexts.empty())
    {
        // Add CPU context
        s_testHardwareContexts.emplace_back();

        // Add GPU contexts for each detected device
        const auto& testDeviceList = getTestHardwareEnvironment()->getTestDeviceList();
        for (const auto& testDevice : testDeviceList)
        {
            s_testHardwareContexts.emplace_back(testDevice.get());
        }
    }

    return s_testHardwareContexts;
}

std::vector<const TestHardwareContext*> getHardwareContextsWithCapability(bool hardwareHasCapability)
{
    std::vector<const TestHardwareContext*> contexts;
    for (const auto& context : getTestHardwareContexts())
    {
        // Include CPU context always, GPU contexts only if capability is supported
        if (!context.isGpuTest() || hardwareHasCapability)
        {
            contexts.push_back(&context);
        }
    }
    return contexts;
}


} // namespace test
} // namespace gmx
