/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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

#include "gromacs/utility/basedefinitions.h"

struct gmx_gpu_info_t;
struct gmx_gpu_opt_t;
struct t_commrec;

namespace gmx
{
class MDLogger;
}

/*! \brief
 * Assigns GPUs to this rank (or checks the manually selected GPU IDs); fills in gpu_opt->dev_use array;
 * initializes all the GPUs of this rank.
 */
void gmx_select_rank_gpu_ids(const gmx::MDLogger &mdlog, const t_commrec *cr,
                             const gmx_gpu_info_t *gpu_info,
                             bool forceUseGPU, bool useGpuNB, bool useGpuPME,
                             gmx_gpu_opt_t *gpu_opt);

/*! \brief
 * Assigns GPUs of this rank to the tasks.
 * Should be called after gmx_select_rank_gpu_ids.
 */
void gmx_select_tasks_gpu_ids(const t_commrec     *cr,
                              gmx_gpu_opt_t       *gpu_opt,
                              bool                 useGpuNB,
                              bool                 useGpuPME);

//! Tasks that can be run on a GPU
enum class GpuTask
{
    NB, PME
};

#define INVALID_GPU_INDEX -1

/*! \internal \brief
 *  Storage of a rank-local information on GPUs assignment to tasks. This stores indices into gpu_opt->dev_use.
 */
class GpuTaskManager
{
    private:
        std::map<GpuTask, int> GpuIndicesByTasks_; /* The assumption is max 1 GPU per task on a single rank. */
    public:
        //! Returns the GPU index into gpu_opt->dev_use for the given GPU task (INVALID_GPU_INDEX if it is not set).
        int gpuIndex(GpuTask task)
        {
            return GpuIndicesByTasks_.count(task) ? GpuIndicesByTasks_[task] : INVALID_GPU_INDEX;
        }
        //! Sets the GPU index into gpu_opt->dev_use for the given GPU task.
        void setGpuIndex(GpuTask task, int gpuIndex)
        {
            GpuIndicesByTasks_[task] = gpuIndex;
        }
};

#endif

