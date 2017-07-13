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
 * Prints GPU information strings on this node into the stderr and log.
 * Only used for logging errors in heterogenous MPI configurations.
 */
static void print_gpu_detection_stats(const gmx::MDLogger  &mdlog,
                                      const gmx_gpu_info_t &gpu_info)
{
    char onhost[HOSTNAMELEN+10];
    int  ngpu;

    if (!gpu_info.bDetectGPUs)
    {
        /* We skipped the detection, so don't print detection stats */
        return;
    }

    ngpu = gpu_info.n_dev;

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

/*! \brief Return whether all selected GPUs are compatible.
 *
 * Given the list of selected GPU device IDs in \c gpu_opt and
 * detected GPUs in \c gpu_info, return whether all selected GPUs are
 * compatible. If not, place a suitable string in \c errorMessage.
 *
 * \param[in]   gpu_info      Information about detected GPUs
 * \param[in]   gpu_opt       pointer to structure holding GPU options
 * \param[out]  errorMessage  pointer to string to hold a possible error message (is not updated when returning true)
 * \returns                   true if every requested GPU is compatible
 */
static bool checkGpuSelection(const gmx_gpu_info_t &gpu_info,
                              const gmx_gpu_opt_t  *gpu_opt,
                              std::string          *errorMessage)
{
    bool        allOK   = true;
    std::string message = "Some of the requested GPUs do not exist, behave strangely, or are not compatible:\n";
    for (int i = 0; i < gpu_opt->n_dev_use; i++)
    {
        GMX_ASSERT(gpu_opt, "Invalid gpu_opt");
        GMX_ASSERT(gpu_opt->dev_use, "Invalid gpu_opt->dev_use");

        int id     = gpu_opt->dev_use[i];
        if (!isGpuCompatible(gpu_info, id))
        {
            allOK    = false;
            message += gmx::formatString("    GPU #%d: %s\n",
                                         id,
                                         getGpuCompatibilityDescription(gpu_info, id));
        }
    }
    if (!allOK && errorMessage)
    {
        *errorMessage = message;
    }
    return allOK;
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

void gmx_select_rank_gpu_ids(const gmx::MDLogger &mdlog, const t_commrec *cr,
                             const gmx_gpu_info_t &gpu_info,
                             bool userSetGpuIds,
                             gmx_gpu_opt_t *gpu_opt)
{
    if (!(cr->duty & DUTY_PP))
    {
        /* Our rank is not doing PP, we don't use a GPU */
        return;
    }

    if (userSetGpuIds)
    {
        /* Check the GPU IDs passed by the user.
         * (GPU IDs have been parsed by gmx_parse_gpu_ids before)
         */
        std::string errorMessage;
        if (!checkGpuSelection(gpu_info, gpu_opt, &errorMessage))
        {
            const bool canHaveHeterogeneousNodes = GMX_LIB_MPI && PAR(cr);
            if (canHaveHeterogeneousNodes)
            {
                print_gpu_detection_stats(mdlog, gpu_info);
            }
            gmx_fatal(FARGS, errorMessage.c_str());
        }
    }
    else
    {
        auto compatibleGpus = getCompatibleGpus(gpu_info);
        assign_rank_gpu_ids(compatibleGpus, gpu_opt, cr->nrank_pp_intranode, cr->rank_pp_intranode);
    }
}
