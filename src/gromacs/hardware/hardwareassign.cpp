/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "hardwareassign.h"

#include "config.h"

#include <cstring>

#include <algorithm>
#include <string>

#include "thread_mpi/threads.h"

#include "gromacs/gmxlib/network.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

#define HOSTNAMELEN 80

/*! \internal \brief
 * Prints GPU information strings on this node into the stderr and log.
 * Only used for logging errors in heterogenous MPI configurations.
 */
static void print_gpu_detection_stats(const gmx::MDLogger  &mdlog,
                                      const gmx_gpu_info_t *gpu_info)
{
    char onhost[HOSTNAMELEN+10];
    int  ngpu;

    if (!gpu_info->bDetectGPUs)
    {
        /* We skipped the detection, so don't print detection stats */
        return;
    }

    ngpu = gpu_info->n_dev;

    /* We only print the detection on one, of possibly multiple, nodes */
    std::strncpy(onhost, " on host ", 10);
    gmx_gethostname(onhost + 9, HOSTNAMELEN);

    if (ngpu > 0)
    {
        std::string gpuDesc = sprint_gpus(gpu_info);
        GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                "%d GPU%s detected%s:\n%s",
                ngpu, (ngpu > 1) ? "s" : "", onhost, gpuDesc.c_str());
    }
    else
    {
        GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted("No GPUs detected%s", onhost);
    }
    // FIXME: This currently only logs on the master rank, which defeats the purpose.
    // A new MDLogger option is required for printing to stderr on all ranks.
    // There is also a question of MPI reduction of the outputs, see Redmine issue #1505.
}

void GpuTaskAssignmentManager::communicateGpuTasksCounts()
{
    const int dev_use_count = tasksToAssign_.size();
    // default values which are valid for a serial case
    int       dev_use_index      = 0;
    int       dev_use_count_node = dev_use_count;
    if ((PAR(cr_) || MULTISIM(cr_))) // Checking for MPI to be initialized
    {
#if GMX_MPI
        const int nrank_intranode = cr_->nrank_intranode;
        const int rank_intranode  = cr_->rank_intranode;
        MPI_Comm  physicalnode_comm;
        MPI_Comm_split(MPI_COMM_WORLD, gmx_physicalnode_id_hash(), rank_intranode, &physicalnode_comm);
        /* Prefix sums of the per-rank GPU task counts */
        MPI_Scan(const_cast<int *>(&dev_use_count), &dev_use_index, 1, MPI_INT, MPI_SUM, physicalnode_comm);
        /* Getting total amount of GPU tasks on this node - the last prefix sum is full sum */
        dev_use_count_node = dev_use_index;
        MPI_Bcast(&dev_use_count_node, 1, MPI_INT, nrank_intranode - 1, physicalnode_comm);
        /* MPI_Scan is inclusive prefix sum, we need exclusive for the starting indices */
        dev_use_index -= dev_use_count;
        MPI_Comm_free(&physicalnode_comm);
#endif
    }
    devUseIndex_      = dev_use_index;
    devUseCountNode_  = dev_use_count_node;
}

void GpuTaskAssignmentManager::assignRankGpuIds()
{
    const int dev_use_count_node = devUseCountNode_;

    GMX_RELEASE_ASSERT(dev_use_count_node >= 1,
                       gmx::formatString("Invalid limit (%d) for the number of GPUs (detected %d compatible GPUs)",
                                         dev_use_count_node, gpuOpt_->n_dev_compatible).c_str());

    if (gpuOpt_->n_dev_compatible == 0)
    {
        char host[HOSTNAMELEN];

        gmx_gethostname(host, HOSTNAMELEN);
        gmx_fatal(FARGS, "A GPU was requested on host %s, but no compatible GPUs were detected. If you intended to use GPU acceleration in a parallel run, you can either avoid using the nodes that don't have GPUs or place CPU tasks on these nodes.", host);
    }

    int nshare = 1; /* Max. number of GPU tasks sharing a single GPU */
    if (dev_use_count_node > gpuOpt_->n_dev_compatible)
    {
        if (dev_use_count_node % gpuOpt_->n_dev_compatible == 0)
        {
            nshare = gmx_gpu_sharing_supported() ? (dev_use_count_node / gpuOpt_->n_dev_compatible) : 1;
        }
        else
        {
            const bool firstGpuRank = (devUseIndex_ == 0) && !tasksToAssign_.empty();
            if (firstGpuRank) // Printing the error message only on one rank
            {
                gmx_fatal(FARGS, "The number of GPU tasks (%d) is not a multiple of the actual number of GPUs (%d). Select a different number of MPI ranks or use the -gpu_id option to manually specify the GPUs to be used.",
                          dev_use_count_node, gpuOpt_->n_dev_compatible);
            }

#if GMX_MPI
            /* We use a global barrier to prevent ranks from continuing with
             * an invalid setup.
             */
            MPI_Barrier(MPI_COMM_WORLD);
#endif
        }
    }

    /* Here we will waste GPUs when dev_use_count_node < gpu_opt->n_dev_compatible */
    gpuOpt_->n_dev_use = std::min(gpuOpt_->n_dev_compatible * nshare, dev_use_count_node);
    if (!gmx_multiple_gpu_per_node_supported())
    {
        gpuOpt_->n_dev_use = std::min(gpuOpt_->n_dev_use, 1);
    }
    snew(gpuOpt_->dev_use, gpuOpt_->n_dev_use);
    for (int i = 0; i != gpuOpt_->n_dev_use; ++i)
    {
        /* TODO: improve this implementation: either sort GPUs or remove the weakest here */
        gpuOpt_->dev_use[i] = gpuOpt_->dev_compatible[i / nshare];
    }
}

void GpuTaskAssignmentManager::registerGpuTask(GpuTask task)
{
    tasksToAssign_.insert(task);
}

void GpuTaskAssignmentManager::selectRankGpuIds(const gmx::MDLogger  &mdlog,
                                                bool                  forceUseGPU)
{
    /* Thsi should be called after all the registerGpuTask() calls,
     * so we know how many GPUs this process can use at most.
     * The actual used GPU count at any point can potentially be smaller.
     */
    communicateGpuTasksCounts();

    if (tasksToAssign_.empty())
    {
        /* Ignore (potentially) manually selected GPUs */
        gpuOpt_->n_dev_use = 0;
        return;
    }

    int              i;
    char             sbuf[STRLEN], stmp[STRLEN];

    /* Bail if binary is not compiled with GPU acceleration, but this is either
     * explicitly (-nb gpu) or implicitly (gpu ID passed) requested. */
    if (forceUseGPU && (GMX_GPU == GMX_GPU_NONE))
    {
        gmx_fatal(FARGS, "GPU acceleration requested, but %s was compiled without GPU support!",
                  gmx::getProgramContext().displayName());
    }

    if (gpuOpt_->bUserSet)
    {
        /* Check the GPU IDs passed by the user.
         * (GPU IDs have been parsed by gmx_parse_gpu_ids before)
         */
        int *checkres;
        snew(checkres, gpuOpt_->n_dev_use);

        int res = check_selected_gpus(checkres, gpuInfo_, gpuOpt_);

        if (!res)
        {
            const bool canHaveHeterogeneousNodes = GMX_LIB_MPI && PAR(cr_);
            if (canHaveHeterogeneousNodes)
            {
                print_gpu_detection_stats(mdlog, gpuInfo_);
            }

            sprintf(sbuf, "Some of the requested GPUs do not exist, behave strangely, or are not compatible:\n");
            for (i = 0; i < gpuOpt_->n_dev_use; i++)
            {
                if (checkres[i] != egpuCompatible)
                {
                    sprintf(stmp, "    GPU #%d: %s\n",
                            gpuOpt_->dev_use[i],
                            gpu_detect_res_str[checkres[i]]);
                    strcat(sbuf, stmp);
                }
            }
            gmx_fatal(FARGS, "%s", sbuf);
        }

        sfree(checkres);
    }
    else if (getenv("GMX_EMULATE_GPU") == nullptr)
    {
        pick_compatible_gpus(gpuInfo_, gpuOpt_);
        /* Assign GPUs to ranks automatically. Intra-rank GPU to task assignment happens in selectTasksGpuIds. */
        assignRankGpuIds();
    }

    /* If the user asked for a GPU, check whether we have a GPU */
    if (forceUseGPU && gpuInfo_->n_dev_compatible == 0)
    {
        gmx_fatal(FARGS, "GPU acceleration requested, but no compatible GPUs were detected.");
    }
}

void GpuTaskAssignmentManager::selectTasksGpuIds()
{
    /* Here we assign indices into gpu_opt->dev_use to the GPU tasks of this rank.
     * gpu_opt->dev_use had already been filled in selectRankGpuIds.
     */
    GMX_RELEASE_ASSERT(gpuOpt_->gpuTasks == nullptr, "GPU task manager should only be initialized once");
    gpuOpt_->gpuTasks.reset(new GpuTaskManager(gpuInfo_, gpuOpt_));
    int gpuIndex = devUseIndex_;
    if (tasksToAssign_.count(GpuTask::NB) > 0)
    {
        gpuOpt_->gpuTasks->setGpuIndex(GpuTask::NB, gpuIndex);
        gpuIndex++;
    }
}

gmx_device_info_t *GpuTaskManager::gpuInfo(GpuTask task)
{
    const size_t deviceInfoSize = sizeof_gpu_dev_info();
    // currently device ID is the same as the index into gpu_dev
    const int    deviceIndex = gpuId(task);
    //TODO this could just live in the gpu_utils to avoid the casts
    return reinterpret_cast<gmx_device_info_t *>(reinterpret_cast<char *>(gpuInfo_->gpu_dev) + deviceIndex * deviceInfoSize);
}

int GpuTaskManager::gpuId(GpuTask task)
{
    int gpuId;
    try
    {
        gpuId = gpuIdsByTasks_[task];
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    return gpuId;
}

//! Sets the GPU index into gpu_opt->dev_use for the given GPU task.
void GpuTaskManager::setGpuIndex(GpuTask task, int gpuIndex)
{
    //TODO this function could probably  accept gpuId instead?
    if (gpuIndex < 0 || gpuIndex >= gpuOpt_->n_dev_use)
    {
        std::string errorString = gmx::formatString("Trying to assign a non-existent GPU: "
                                                    "there are %d %s-selected GPU(s), but #%d was requested.",
                                                    gpuOpt_->n_dev_use, gpuOpt_->bUserSet ? "user" : "auto", gpuIndex);
        gmx_incons(errorString.c_str());
    }

    const int gpuId    = get_gpu_device_id(gpuInfo_, gpuOpt_, gpuIndex);
    gpuIdsByTasks_[task] = gpuId;

}
