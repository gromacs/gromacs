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

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/hardwareassign.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/loggerbuilder.h"

namespace gmx
{
namespace test
{
//! This constructs the test environment
PmeTestEnvironment * const pmeEnv = (PmeTestEnvironment *)::testing::AddGlobalTestEnvironment(new PmeTestEnvironment);

std::string GpuTestHardwareContext::getDescription()
{
    GMX_RELEASE_ASSERT(get_current_cuda_gpu_device_id() == id_, "The GPU context got unexpectedly changed");
    char stmp[200] = {};
    get_gpu_device_info_string(stmp, &pmeEnv->getHardwareInfo()->gpu_info, id_);
    return std::string(stmp);
}

void GpuTestHardwareContext::activate()
{
#if GMX_GPU != GMX_GPU_CUDA
    GMX_THROW(NotImplementedError("Not implemented"));
#endif
    set_current_cuda_gpu_device_id(id_);
}

//! Simple hardware initialization
void PmeTestEnvironment::hardwareInit()
{
    commrec_.reset(init_commrec());
    gmx_init_intranode_counters(commrec_.get());
    gpuOptions_.reset(new gmx_gpu_opt_t {});
    LoggerBuilder builder;
    LoggerOwner   logOwner(builder.build());
    MDLogger      log(logOwner.logger());
    hardwareInfo_.reset(gmx_detect_hardware(log, commrec_.get(), true));
    pick_compatible_gpus(&hardwareInfo_.get()->gpu_info, gpuOptions_.get());
    gmx_select_rank_gpu_ids(log, commrec_.get(), &hardwareInfo_.get()->gpu_info, false, gpuOptions_.get());
}

void PmeTestEnvironment::SetUp()
{
    auto emptyContext = std::make_shared<EmptyTestHardwareContext>();
    hardwareContextsByMode_[CodePath::CPU]  = {emptyContext};

#if GMX_GPU == GMX_GPU_CUDA
    hardwareInit();  // TODO - move this line out of conditional, when gmx_select_gpu_ids() stops failing in CPU-only reference build
    TestHardwareContexts gpuContexts;
    for (ssize_t i = 0; i < gpuOptions_->n_dev_compatible; i++)
    {
        auto gpuContext = std::make_shared<GpuTestHardwareContext>(GpuTestHardwareContext(gpuOptions_->dev_compatible[i]));
        gpuContexts.push_back(gpuContext);
    }
    hardwareContextsByMode_[CodePath::CUDA] = gpuContexts;
#endif
}

const TestHardwareContexts &GetContextsForMode(CodePath mode)
{
    return pmeEnv->getHardwareContexts(mode);
}

}
}
