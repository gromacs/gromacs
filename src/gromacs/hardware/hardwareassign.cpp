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

#include "gromacs/gmxlib/network.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

#define HOSTNAMELEN 80

/*! \internal \brief
 * This function is responsible for mapping the GPUs to the processes on a single node
 * (filling the gpu_opt->dev_use array).
 *
 * \param[in]        compatibleGpus       Vector of GPUs that are compatible
 * \param[in,out]    gpu_opt              Input/output GPU assignment data.
 * \param[in]        nrank                Number of PP GPU ranks on the node.
 * \param[in]        rank                 Index of PP GPU rank on the node.
 *
 * This selects the GPUs we will use. This is an operation local to each physical node.
 * If we have less MPI ranks than GPUs, we will waste some GPUs.
 */
static void assign_rank_gpu_ids(const std::vector<int> &compatibleGpus,
                                gmx_gpu_opt_t *gpu_opt, int nrank, int rank)
{
    int numCompatibleGpus = static_cast<int>(compatibleGpus.size());
    GMX_RELEASE_ASSERT(gpu_opt, "Invalid gpu_opt pointer passed");
    GMX_RELEASE_ASSERT(nrank >= 1,
                       gmx::formatString("Invalid limit (%d) for the number of GPUs (detected %d compatible GPUs)",
                                         rank, numCompatibleGpus).c_str());

    if (numCompatibleGpus == 0)
    {
        char host[HOSTNAMELEN];

        gmx_gethostname(host, HOSTNAMELEN);
        gmx_fatal(FARGS, "A GPU was requested on host %s, but no compatible GPUs were detected. All nodes with PP ranks need to have GPUs. If you intended to use GPU acceleration in a parallel run, you can either avoid using the nodes that don't have GPUs or place PME ranks on these nodes.", host);
    }

    int nshare;

    nshare = 1;
    if (nrank > numCompatibleGpus)
    {
        if (nrank % numCompatibleGpus == 0)
        {
            nshare = gmx_gpu_sharing_supported() ? nrank/numCompatibleGpus : 1;
        }
        else
        {
            if (rank == 0)
            {
                gmx_fatal(FARGS, "The number of MPI ranks (%d) in a physical node is not a multiple of the number of GPUs (%d). Select a different number of MPI ranks or use the -gpu_id option to manually specify the GPU to be used.",
                          nrank, numCompatibleGpus);
            }

#if GMX_MPI
            /* We use a global barrier to prevent ranks from continuing with
             * an invalid setup.
             */
            MPI_Barrier(MPI_COMM_WORLD);
#endif
        }
    }

    /* Here we will waste GPUs when nrank < numCompatibleGpus */
    gpu_opt->n_dev_use = std::min(numCompatibleGpus*nshare, nrank);
    if (!gmx_multiple_gpu_per_node_supported())
    {
        gpu_opt->n_dev_use = std::min(gpu_opt->n_dev_use, 1);
    }
    snew(gpu_opt->dev_use, gpu_opt->n_dev_use);
    for (int i = 0; i != gpu_opt->n_dev_use; ++i)
    {
        /* TODO: improve this implementation: either sort GPUs or remove the weakest here */
        gpu_opt->dev_use[i] = compatibleGpus[i/nshare];
    }
}

/*! \brief Check that all user-selected GPUs are compatible.
 *
 * Given the list of selected GPU device IDs in \c gpu_opt and
 * detected GPUs in \c gpu_info, gives a fatal error unless all
 * selected GPUs are compatible
 *
 * The error is given with a suitable descriptive message, which will
 * have context if this check is done after the hardware detection
 * results have been reported to the user. However, note that only the
 * GPUs detected on the master rank are reported, because of the
 * existing limitations of that reporting.
 *
 * \todo Note that the selected GPUs can be different on each rank,
 * and the IDs of compatible GPUs can be different on each node, so
 * this routine ought to do communication to determine whether all
 * ranks are able to proceed. Currently this relies on the MPI runtime
 * to kill the other processes because GROMACS lacks the appropriate
 * infrastructure to do a good job of coordinating error messages and
 * behaviour across MPMD ranks and multiple simulations.
 *
 * \param[in]   gpu_info       GPU information including result of compatibility check.
 * \param[in]   gpu_opt        Mapping of GPU IDs derived from the user, e.g. via mdrun -gpu_id.
 */
static void exitUnlessGpuSelectionIsValid(const gmx_gpu_info_t &gpu_info,
                                          const gmx_gpu_opt_t  &gpu_opt)
{
    int         numIncompatibleGpuIds = 0;
    std::string message
        = "Some of the requested GPUs do not exist, behave strangely, or are not compatible:\n";

    for (int index = 0; index < gpu_opt.n_dev_use; ++index)
    {
        int gpuId = gpu_opt.dev_use[index];
        if (!isGpuCompatible(gpu_info, gpuId))
        {
            numIncompatibleGpuIds++;
            message += gmx::formatString("    GPU #%d: %s\n",
                                         gpuId,
                                         getGpuCompatibilityDescription(gpu_info, gpuId));
        }
    }

    if (numIncompatibleGpuIds > 0)
    {
        gmx_fatal(FARGS, message.c_str());
    }
}

std::vector<int> getCompatibleGpus(const gmx_gpu_info_t &gpu_info)
{
    // Possible minor over-allocation here, but not important for anything
    std::vector<int> compatibleGpus;
    compatibleGpus.reserve(gpu_info.n_dev);
    for (int i = 0; i < gpu_info.n_dev; i++)
    {
        GMX_ASSERT(gpu_info.gpu_dev, "Invalid gpu_info.gpu_dev");
        if (isGpuCompatible(gpu_info, i))
        {
            compatibleGpus.push_back(i);
        }
    }
    return compatibleGpus;
}

void mapPpRanksToGpus(bool                  rankCanUseGpu,
                      const t_commrec      *cr,
                      const gmx_gpu_info_t &gpu_info,
                      bool                  userSetGpuIds,
                      gmx_gpu_opt_t        *gpu_opt)
{
    if (!rankCanUseGpu)
    {
        return;
    }

    if (userSetGpuIds)
    {
        exitUnlessGpuSelectionIsValid(gpu_info, *gpu_opt);
    }
    else
    {
        auto compatibleGpus = getCompatibleGpus(gpu_info);
        assign_rank_gpu_ids(compatibleGpus, gpu_opt, cr->nrank_pp_intranode, cr->rank_pp_intranode);
    }
}
