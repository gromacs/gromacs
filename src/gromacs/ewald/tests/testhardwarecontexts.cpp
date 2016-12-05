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
#include "gromacs/taskassignment/hardwareassign.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/loggerbuilder.h"

namespace gmx
{
namespace test
{

//! This constructs the test environment
PmeTestEnvironment * const pmeEnv = (PmeTestEnvironment *)::testing::AddGlobalTestEnvironment(new PmeTestEnvironment);

//! Simple hardware initialization
void PmeTestEnvironment::hardwareInit()
{
    unique_cptr<t_commrec, done_commrec> commrec(init_commrec());
    gmx_init_intranode_counters(commrec.get());
    LoggerBuilder builder;
    LoggerOwner   logOwner(builder.build());
    MDLogger      log(logOwner.logger());
    hardwareInfo_.reset(gmx_detect_hardware(log, commrec.get()));
}

void PmeTestEnvironment::SetUp()
{
    TestHardwareContext emptyContext("", nullptr);
    hardwareContextsByMode_[CodePath::CPU].push_back(emptyContext);

    hardwareInit();

    // Constructing contexts for all compatible GPUs - will be empty on non-GPU builds
    TestHardwareContexts gpuContexts;
    const auto           compatibleGpus = getCompatibleGpus(hardwareInfo_->gpu_info);
    for (int gpuIndex : compatibleGpus)
    {
        char        stmp[200] = {};
        get_gpu_device_info_string(stmp, hardwareInfo_->gpu_info, gpuIndex);
        std::string description = "(GPU " + std::string(stmp) + ") ";
        auto       *gpuInfo     = reinterpret_cast<gmx_device_info_t *>(reinterpret_cast<char *>(hardwareInfo_->gpu_info.gpu_dev) + gpuIndex * sizeof_gpu_dev_info());
        //TODO move previous line to gpu_utils and reuse in hardwareassign.cpp
        gpuContexts.emplace_back(TestHardwareContext(description.c_str(), gpuInfo));
    }
#if GMX_GPU == GMX_GPU_CUDA
    hardwareContextsByMode_[CodePath::CUDA] = gpuContexts;
#endif
}

}
}
