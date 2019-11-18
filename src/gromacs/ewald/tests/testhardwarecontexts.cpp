/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019, by the GROMACS development team, led by
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
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include "testhardwarecontexts.h"

#include <memory>

#include "gromacs/ewald/pme.h"
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

TestHardwareContext::~TestHardwareContext() = default;

const char* codePathToString(CodePath codePath)
{
    switch (codePath)
    {
        case CodePath::CPU: return "CPU";
        case CodePath::GPU: return "GPU";
        default: GMX_THROW(NotImplementedError("This CodePath should support codePathToString"));
    }
}

/* Implements the "construct on first use" idiom to avoid any static
 * initialization order fiasco.
 *
 * Note that thread-safety of the initialization is guaranteed by the
 * C++11 language standard.
 *
 * The pointer itself (not the memory it points to) has no destructor,
 * so there is no deinitialization issue.  See
 * https://isocpp.org/wiki/faq/ctors for discussion of alternatives
 * and trade-offs. */
const PmeTestEnvironment* getPmeTestEnv()
{
    static PmeTestEnvironment* pmeTestEnvironment = nullptr;
    if (pmeTestEnvironment == nullptr)
    {
        // Ownership of the TestEnvironment is taken by GoogleTest, so nothing can leak
        pmeTestEnvironment = static_cast<PmeTestEnvironment*>(
                ::testing::AddGlobalTestEnvironment(new PmeTestEnvironment));
    }
    return pmeTestEnvironment;
}

void callAddGlobalTestEnvironment()
{
    getPmeTestEnv();
}

//! Simple hardware initialization
static gmx_hw_info_t* hardwareInit()
{
    PhysicalNodeCommunicator physicalNodeComm(MPI_COMM_WORLD, gmx_physicalnode_id_hash());
    return gmx_detect_hardware(MDLogger{}, physicalNodeComm);
}

void PmeTestEnvironment::SetUp()
{
    hardwareContexts_.emplace_back(std::make_unique<TestHardwareContext>(CodePath::CPU, "", nullptr));

    hardwareInfo_ = hardwareInit();
    if (!pme_gpu_supports_build(nullptr) || !pme_gpu_supports_hardware(*hardwareInfo_, nullptr))
    {
        // PME can only run on the CPU, so don't make any more test contexts.
        return;
    }
    // Constructing contexts for all compatible GPUs - will be empty on non-GPU builds
    for (int gpuIndex : getCompatibleGpus(hardwareInfo_->gpu_info))
    {
        const gmx_device_info_t* deviceInfo = getDeviceInfo(hardwareInfo_->gpu_info, gpuIndex);
        init_gpu(deviceInfo);

        char stmp[200] = {};
        get_gpu_device_info_string(stmp, hardwareInfo_->gpu_info, gpuIndex);
        std::string description = "(GPU " + std::string(stmp) + ") ";
        hardwareContexts_.emplace_back(std::make_unique<TestHardwareContext>(
                CodePath::GPU, description.c_str(), deviceInfo));
    }
}

} // namespace test
} // namespace gmx
