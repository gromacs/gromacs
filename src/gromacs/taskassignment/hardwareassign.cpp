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
#include <set>
#include <string>
#include <vector>

#include "gromacs/gmxlib/network.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

#define HOSTNAMELEN 80

namespace gmx
{

std::vector<int> parseGpuTaskAssignment(const std::string &gpuTaskAssignment)
{
    std::vector<int> digits;
    if (gpuTaskAssignment.empty())
    {
        return digits;
    }

    /* Parse a "plain" or comma-separated GPU ID string which contains
     * a sequence of digits corresponding to GPU IDs; the order will
     * indicate the assignment of GPU tasks on this node to GPU
     * device IDs on this node. */
    try
    {
        digits = parseDigitsFromString(gpuTaskAssignment);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    if (digits.empty())
    {
        gmx_fatal(FARGS, "Empty GPU ID string encountered.\n"
                  "An empty, delimiter-free, or comma-separated sequence of valid numeric IDs of available GPUs is required.\n");
    }
    return digits;
}

/*! \brief This function is responsible for the automated mapping the
 * GPUs to the processes on a single node.
 *
 * This selects the GPUs we will use. This is an operation local to each physical node.
 * If we have less MPI ranks than GPUs, we will waste some GPUs.
 *
 * \param[in]        compatibleGpus       Vector of GPUs that are compatible
 * \param[in]        nrank                Number of PP GPU ranks on the node.
 * \param[in]        rank                 Index of PP GPU rank on the node.
 *
 * \returns The assignment of GPU tasks on ranks of this node to GPU devices on this node.
 */
static std::vector<int> assign_rank_gpu_ids(const std::vector<int> &compatibleGpus,
                                            int nrank, int rank)
{
    int numCompatibleGpus = static_cast<int>(compatibleGpus.size());
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
            nshare = nrank/numCompatibleGpus;
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
    std::vector<int> taskAssignment;
    taskAssignment.resize(std::min(numCompatibleGpus*nshare, nrank));
    for (size_t i = 0; i != taskAssignment.size(); ++i)
    {
        /* TODO: improve this implementation: either sort GPUs or remove the weakest here */
        taskAssignment[i] = compatibleGpus[i/nshare];
    }
    return taskAssignment;
}

/*! \brief Check that all user-selected GPUs are compatible.
 *
 * Given the \c userGpuTaskAssignment and \c compatibleGPUs, give a fatal
 * error if any selected GPUs is not compatible
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
 * \param[in]   gpu_info               GPU information including device description.
 * \param[in]   compatibleGpus         Vector of compatible GPUs
 * \param[in]   userGpuTaskAssignment  The GPU selection from the user.
 */
static void exitUnlessUserGpuTaskAssignmentIsValid(const gmx_gpu_info_t   &gpu_info,
                                                   const std::vector<int> &compatibleGpus,
                                                   const std::vector<int> &userGpuTaskAssignment)
{
    int         numIncompatibleGpuIds = 0;
    std::string message
        = "Some of the requested GPUs do not exist, behave strangely, or are not compatible:\n";

    for (const auto &gpuId : userGpuTaskAssignment)
    {
        if (std::find(compatibleGpus.begin(), compatibleGpus.end(), gpuId) == compatibleGpus.end())
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

std::vector<int> mapPpRanksToGpus(bool                    rankCanUseGpu,
                                  const t_commrec        *cr,
                                  const gmx_gpu_info_t   &gpu_info,
                                  const gmx_hw_opt_t     &hw_opt)
{
    std::vector<int> taskAssignment;

    if (!rankCanUseGpu)
    {
        return taskAssignment;
    }

    auto compatibleGpus = getCompatibleGpus(gpu_info);
    if (!hw_opt.gpuIdTaskAssignment.empty())
    {
        auto userGpuTaskAssignment = parseGpuTaskAssignment(hw_opt.gpuIdTaskAssignment);
        exitUnlessUserGpuTaskAssignmentIsValid(gpu_info, compatibleGpus, userGpuTaskAssignment);
        taskAssignment = userGpuTaskAssignment;
    }
    else
    {
        taskAssignment = assign_rank_gpu_ids(compatibleGpus, cr->nrank_pp_intranode, cr->rank_pp_intranode);
    }
    return taskAssignment;
}

} // namespace

/*! \brief Return the number of PP rank pairs that share a GPU device between them.
 *
 * Sharing GPUs among multiple PP ranks is possible via either user or
 * automated selection. */
static int gmx_count_gpu_dev_shared(const std::vector<int> &gpuTaskAssignment,
                                    bool                    userSetGpuIds)
{
    int      same_count    = 0;

    if (userSetGpuIds)
    {
        GMX_RELEASE_ASSERT(!gpuTaskAssignment.empty(),
                           "The user cannot choose an empty set of GPU IDs, code is wrong somewhere");
        size_t ngpu = gpuTaskAssignment.size();

        for (size_t i = 0; i < ngpu - 1; i++)
        {
            for (size_t j = i + 1; j < ngpu; j++)
            {
                same_count      += (gpuTaskAssignment[i] ==
                                    gpuTaskAssignment[j]);
            }
        }
    }

    return same_count;
}

/* Count and return the number of unique GPUs (per node) selected.
 *
 * As sharing GPUs among multiple PP ranks is possible, the number of
 * GPUs used (per node) can be different from the number of GPU IDs
 * used.
 */
static size_t gmx_count_gpu_dev_unique(const std::vector<int> &gpuTaskAssignment)
{
    std::set<int> uniqIds;
    for (const auto &deviceId : gpuTaskAssignment)
    {
        uniqIds.insert(deviceId);
    }
    return uniqIds.size();
}

void reportGpuUsage(const gmx::MDLogger    &mdlog,
                    const gmx_gpu_info_t   &gpu_info,
                    bool                    userSetGpuIds,
                    const std::vector<int> &gpuTaskAssignment,
                    size_t                  numPpRanks,
                    bool                    bPrintHostName)
{
    if (gpuTaskAssignment.empty())
    {
        return;
    }

    std::string output;
    {
        std::string gpuIdsString =
            formatAndJoin(gpuTaskAssignment, ",", gmx::StringFormatter("%d"));
        size_t      numGpusInUse = gmx_count_gpu_dev_unique(gpuTaskAssignment);
        bool        bPluralGpus  = numGpusInUse > 1;

        if (bPrintHostName)
        {
            char host[STRLEN];
            gmx_gethostname(host, STRLEN);
            output += gmx::formatString("On host %s ", host);
        }
        output += gmx::formatString("%zu GPU%s %sselected for this run.\n"
                                    "Mapping of GPU ID%s to the %d PP rank%s in this node: %s\n",
                                    numGpusInUse, bPluralGpus ? "s" : "",
                                    userSetGpuIds ? "user-" : "auto-",
                                    bPluralGpus ? "s" : "",
                                    numPpRanks,
                                    (numPpRanks > 1) ? "s" : "",
                                    gpuIdsString.c_str());
    }

    int same_count = gmx_count_gpu_dev_shared(gpuTaskAssignment, userSetGpuIds);

    if (same_count > 0)
    {
        output += gmx::formatString("NOTE: You assigned %s to multiple ranks.\n",
                                    same_count > 1 ? "GPU IDs" : "a GPU ID");
    }

    if (static_cast<size_t>(gpu_info.n_dev_compatible) > numPpRanks)
    {
        /* TODO In principle, this warning could be warranted only on
         * ranks on some nodes, but we lack the infrastructure to do a
         * good job of reporting that. */
        output += gmx::formatString("NOTE: potentially sub-optimal launch configuration using fewer\n"
                                    "      PP ranks on a node than GPUs available on that node.\n");
    }

    /* NOTE: this print is only for and on one physical node */
    GMX_LOG(mdlog.warning).appendText(output);
}
