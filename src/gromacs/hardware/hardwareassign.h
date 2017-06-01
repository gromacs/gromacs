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
 * A handle/description of a GPU context.
 * \todo Move routines such as get_gpu_device_info_string() into this.
 */
struct GpuContext
{
    //! A platform-dependent device information pointer
    gmx_device_info_t *gpuInfo_;
    //! A numerical identificator of the GPU - expected to be used for human-readable output mostly
    int                gpuId_;
};

//! Assignment of GPUs to tasks on the rank. Assumption is max 1 GPU context per task type.
typedef std::map<GpuTask, GpuContext> GpuContextsMap;

/*! \libinternal \brief
 * Handles setup-time assignment of GPUs to tasks.
 * To use this class one has to call registerGpuTask() for each type of task desired;
 * then selectRankGpus() and finally selectTasksGpus()
 * (the last call returns the mapping of GPU context handles to task types).
 * This is all done in the convenience wrapper createGpuAssignment(), which is declared below.
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

    /*! \brief
     * Creates and returns the GPU context description, based on index into gpu_opt->dev_use.
     * Will terminate the program if \p gpuIndex is out of bounds.
     */
    GpuContext getGpuContext(int gpuIndex);

    public:
        //! Default constructor is disallowed as the GPU and communication structures are required
        GpuTaskAssignmentManager() = delete;
        //! Constructs the assignment manager
        GpuTaskAssignmentManager(const gmx_gpu_info_t *gpuInfo, gmx_gpu_opt_t *gpuOpt) :
            devUseIndex_(0), devUseCountNode_(0), gpuInfo_(gpuInfo), gpuOpt_(gpuOpt){}
        /*! \brief
         * Registers the intent of the task to be run on the GPU.
         * Will terminate the program on the non-GPU build.
         */
        void registerGpuTask(GpuTask task);
        /*! \brief
         * Assigns GPUs to this rank (or checks the manually selected GPU IDs); fills in gpu_opt->dev_use array.
         * Uses cr to communicate the GPU taks counts intra-node.
         */
        void selectRankGpus(const gmx::MDLogger &mdlog, const t_commrec *cr);
        /*! \brief
         * Assigns GPUs of this rank to the tasks and returns the resulting mapping. Should be called after selectRankGpus().
         */
        GpuContextsMap selectTasksGpus();
};

/*! \libinternal \brief
 * Storage of a rank-local information on GPUs assignment to tasks.
 * Provides a convenience wrapper over GpuContextsMap.
 */
class GpuTaskManager
{
    private:
        //! Underlying storage.
        GpuContextsMap gpuContextsByTasks_;
    public:
        //! Default constructor is disallowed
        GpuTaskManager() = delete;
        //! Constructs the object from the underlying map
        GpuTaskManager(const GpuContextsMap &gpuContextsByTasks) : gpuContextsByTasks_(gpuContextsByTasks){}
        //! Returns the GPU ID for the given GPU task
        //! \throws std::out_of_range if there is no such GPU task assigned on this rank
        //! \returns GPU ID of the task if it was assigned
        int gpuId(GpuTask task) const;
        //! Returns the pointer to the platform-dependent GPU info for the given GPU task
        //! \returns The pointer to the GPU info; nullptr if such GpuTask was not registered on this rank.
        gmx_device_info_t *gpuInfo(GpuTask task) const;
        //! Returns the number of GPU tasks on this rank
        size_t rankGpuTasksCount() const;
        //! Returns all the GPU contexts used by this rank
        std::set<gmx_device_info_t *> rankGpuInfos() const;
};

//! A wrapper function which assigns GPUs to tasks using the classes above
GpuTaskManager createGpuAssignment(const gmx::MDLogger &mdlog, const t_commrec *cr,
                                   const gmx_gpu_info_t &gpuInfo, gmx_gpu_opt_t &gpuOpt,
                                   bool useGpuNB, bool useGpuPME);

#endif
