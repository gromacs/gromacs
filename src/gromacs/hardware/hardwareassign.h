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
    NB,
    PME
};

/*! \libinternal \brief Select the compatible GPUs
 *
 * This function filters gpu_info->gpu_dev for compatible gpus based
 * on the previously run compatibility tests. Sets
 * gpu_info->dev_compatible and gpu_info->n_dev_compatible.
 *
 * \param[in]     gpu_info    pointer to structure holding GPU information
 * \param[out]    gpu_opt     pointer to structure holding GPU options
 */
void pickCompatibleGpus(const gmx_gpu_info_t *gpu_info,
                        gmx_gpu_opt_t        *gpu_opt);

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
    //! Intermediate GPU assignment information that is partially modified by this object
    gmx_gpu_opt_t        *gpuOpt_;
    //! GPU mapping to fill in
    class GpuTaskManager *gpuTasks_;
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
        GpuTaskAssignmentManager(const gmx_gpu_info_t *gpuInfo, gmx_gpu_opt_t *gpuOpt, GpuTaskManager *gpuTasks) :
            devUseIndex_(0), devUseCountNode_(0), gpuInfo_(gpuInfo), gpuOpt_(gpuOpt), gpuTasks_(gpuTasks){}

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
        //! TODO: Consider moving this to gpu_utils/removing
        //! \throws std::out_of_range if there is no such GPU task assigned on this rank
        //! \returns GPU ID of the task if it was assigned
        int gpuId(GpuTask task) const;
        //! Returns the pointer to the platform-dependent GPU info for the given GPU task
        gmx_device_info_t *gpuInfo(GpuTask task) const;
        //! Sets the GPU index into gpu_opt->dev_use for the given GPU task.
        void setGpuIndex(GpuTask task, int gpuIndex);
        //! Tells if this rank uses any GPUs
        bool rankHasGpuTasks() const;
};

#endif
