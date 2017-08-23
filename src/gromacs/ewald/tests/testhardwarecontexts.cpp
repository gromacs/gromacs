/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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

#include "config.h"

#include "gromacs/gmxlib/network.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/loggerbuilder.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{
namespace test
{

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
const PmeTestEnvironment *getPmeTestEnv()
{
    static PmeTestEnvironment *pmeTestEnvironment = nullptr;
    if (pmeTestEnvironment == nullptr)
    {
        // Ownership of the TestEnvironment is taken by GoogleTest, so nothing can leak
        pmeTestEnvironment = static_cast<PmeTestEnvironment *>(::testing::AddGlobalTestEnvironment(new PmeTestEnvironment));
    }
    return pmeTestEnvironment;
}

void callAddGlobalTestEnvironment()
{
    getPmeTestEnv();
}

//! Simple hardware initialization
static gmx_hw_info_t *hardwareInit()
{
    unique_cptr<t_commrec, done_commrec> commrec(init_commrec());
    gmx_init_intranode_counters(commrec.get());
    LoggerBuilder builder;
    LoggerOwner   logOwner(builder.build());
    MDLogger      log(logOwner.logger());
    return gmx_detect_hardware(log, commrec.get());
}

void PmeTestEnvironment::SetUp()
{
    TestHardwareContext emptyContext("", nullptr);
    hardwareContextsByMode_[CodePath::CPU].push_back(emptyContext);

    hardwareInfo_ = hardwareInit();

    // Constructing contexts for all compatible GPUs - will be empty on non-GPU builds
    TestHardwareContexts gpuContexts;
    for (int gpuIndex : getCompatibleGpus(hardwareInfo_->gpu_info))
    {
        char        stmp[200] = {};
        get_gpu_device_info_string(stmp, hardwareInfo_->gpu_info, gpuIndex);
        std::string description = "(GPU " + std::string(stmp) + ") ";
        gpuContexts.emplace_back(TestHardwareContext(description.c_str(), getDeviceInfo(hardwareInfo_->gpu_info, gpuIndex)));
    }
#if GMX_GPU == GMX_GPU_CUDA
    hardwareContextsByMode_[CodePath::CUDA] = gpuContexts;
#endif
}

void PmeTestEnvironment::TearDown()
{
    gmx_hardware_info_free();
}

}
}
