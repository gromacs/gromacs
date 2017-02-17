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
#ifndef GMX_HARDWARE_HARDWAREASSIGN_H
#define GMX_HARDWARE_HARDWAREASSIGN_H

#include <map>
#include <set>

#include "gromacs/utility/basedefinitions.h"

struct gmx_device_info_t;
struct gmx_gpu_info_t;
struct gmx_gpu_opt_t;
struct t_commrec;

namespace gmx
{
class MDLogger;
}

//! Tasks that can be run on a GPU
enum class GpuTask
{
    NB
};

/*! \libinternal \brief
 * Handles setup-time assignment of GPUs to tasks
 */
class GpuTaskAssignmentManager
{
    //! A starting index into gpu_opt->dev_use GPU ID array for this rank
    int                   devUseIndex_;
    //! Number of GPU tasks on the physical node.
    int                   devUseCountNode_;
    //! GPU information
    const gmx_gpu_info_t *gpuInfo_;
    //! GPU assignment information that is modified by this object
    gmx_gpu_opt_t        *gpuOpt_;
    //! Tasks to assign during selectRank/TasksGpuIds
    std::set<GpuTask>     tasksToAssign_;

    /*! \brief
     * This function is responsible for default mapping of the GPUs to the processes on a single node
     * (filling the gpu_opt->dev_use array).
     *
     * This selects the GPUs we will use. This is an operation local to each physical node.
     * If we have less MPI ranks than GPUs, we will waste some GPUs.
     */
    void assignRankGpuIds();

    public:
        //! Default constructor is disallowed as the GPU and communication structures are required
        GpuTaskAssignmentManager() = delete;
        //! Constructs the assignment manager
        GpuTaskAssignmentManager(const gmx_gpu_info_t *gpuInfo, gmx_gpu_opt_t *gpuOpt) :
            devUseIndex_(0), devUseCountNode_(0), gpuInfo_(gpuInfo), gpuOpt_(gpuOpt){}

        /*! \brief
         * Registers the intent of the task to be run on the GPU
         */
        void registerGpuTask(GpuTask task);

        /*! \brief
         * Assigns GPUs to this rank (or checks the manually selected GPU IDs); fills in gpu_opt->dev_use array.
         * Uses cr to communicate the GPU taks counts intra-node.
         */
        void selectRankGpuIds(const gmx::MDLogger &mdlog, const t_commrec *cr);

        /*! \brief
         * Assigns GPUs of this rank to the tasks. Should be called after selectRankGpuIds.
         */
        void selectTasksGpuIds();
};

/*! \libinternal \brief
 *  Storage of a rank-local information on GPUs assignment to tasks.
 */
class GpuTaskManager
{
    private:
        //! Storing max 1 GPU ID per task type on any rank. */
        std::map<GpuTask, int> gpuIdsByTasks_;
        //! GPU information
        const gmx_gpu_info_t  *gpuInfo_;
        //! GPU assignment information
        const gmx_gpu_opt_t   *gpuOpt_;
    public:
        //! Default constructor is disallowed as the GPU info is required
        GpuTaskManager() = delete;
        //! Constructs the empty mapping
        GpuTaskManager(const gmx_gpu_info_t *gpuInfo, const gmx_gpu_opt_t *gpuOpt) : gpuInfo_(gpuInfo), gpuOpt_(gpuOpt){}
        //! Returns the GPU ID for the given GPU task
        int gpuId(GpuTask task);
        //! Returns the pointer to the platform-dependent GPU info for the given GPU task
        gmx_device_info_t *gpuInfo(GpuTask task);
        //! Sets the GPU index into gpu_opt->dev_use for the given GPU task.
        void setGpuIndex(GpuTask task, int gpuIndex);
};

#endif
