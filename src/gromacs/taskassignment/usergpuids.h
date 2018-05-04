/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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
/*! \defgroup module_taskassignment Assigning simulation tasks to hardware (taskassignment)
 * \ingroup group_mdrun
 * \brief Provides code that manages assignment of simulation tasks to hardware.
 */
/*! \libinternal
 * \file
 * \brief Declares routines for handling user-specified GPU IDs.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_taskassignment
 * \inlibraryapi
 */
#ifndef GMX_TASKASSIGNMENT_USERGPUIDS_H
#define GMX_TASKASSIGNMENT_USERGPUIDS_H

#include <cstddef>

#include <string>
#include <vector>

#include "gromacs/utility/arrayref.h"

struct gmx_gpu_info_t;

namespace gmx
{

class MDLogger;

/*! \brief Parse a GPU ID string into a container describing the task types and associated device IDs.
 *
 * \param[in]   gpuIdString  String like "013" or "0,1,3" typically
 *                           supplied by the user to mdrun -gpu_id or -gputasks.
 *                           Must contain only decimal digits, or only decimal
 *                           digits separated by comma delimiters. A terminal
 *                           comma is accceptable (and required to specify a
 *                           single ID that is larger than 9).
 *
 * \returns  A vector of GPU ID task mappings, like { 0, 1, 3 }
 *
 * \throws   std::bad_alloc     If out of memory.
 *           InvalidInputError  If an invalid character is found (ie not a digit or ',').
 */
std::vector<int>
parseUserGpuIds(const std::string &gpuIdString);

/*! \brief Make a vector containing \c numGpuTasks IDs of the IDs found in \c compatibleGpus.
 *
 * \throws  std::bad_alloc          If out of memory
 *
 * \returns A sorted vector of IDs of compatible vectors, whose
 * length matches that of the number of GPU tasks required.
 */
std::vector<int>
makeGpuIds(ArrayRef<const int> compatibleGpus,
           size_t              numGpuTasks);

/*! \brief Convert a container of GPU deviced IDs to a string that
 * can be used by gmx tune_pme as input to mdrun -gputasks.
 *
 * Produce a valid input for mdrun -gputasks that refers to the device
 * IDs in \c gpuIds but produces a mapping for \c
 * totalNumberOfTasks tasks. Note that gmx tune_pme does not
 * currently support filling mdrun -gputasks.
 *
 * \param[in]   gpuIds              Container of device IDs
 * \param[in]   totalNumberOfTasks  Total number of tasks for the output mapping produced by the returned string.
 *
 * \returns  A string that is suitable to pass to mdrun -gputasks.
 *
 * \throws   std::bad_alloc     If out of memory.
 */
std::string
makeGpuIdString(const std::vector<int> &gpuIds, int totalNumberOfTasks);


/*! \brief Return the set of compatible GPUs that can be used.
 *
 * If the user chose GPU IDs, throw if the GPU IDs could not be
 * parsed, were out of range, duplicated any GPU IDs, or referred to
 * any devices that have unsuitable status. If valid, return the
 * user-selected GPU IDs, otherwise, return all the IDs of available
 * compatible GPUs.
 *
 * The exception thrown has a suitable descriptive message, which will
 * have context if this check is done after the hardware detection
 * results have been reported to the user. However, note that only the
 * GPUs detected on the master rank are reported, because of the
 * existing limitations of that reporting.
 *
 * \param[in]   mdlog                   Logging object to write to.
 * \param[in]   compatibleGpuIds        The detected compatible availble GPU IDs.
 * \param[in]   numGpuDevicesDetected   The total number of GPU devices detected.
 * \param[in]   userGpuIdChoicesString  The GPU IDs selected by the user.
 *
 * \throws  std::bad_alloc     If out of memory.
 *          InvalidInputError  If the user's GPU ID choices are invalid.
 * \returns                    The set of compatible GPU IDs that can be used.
 */
std::vector<int>
determineGpuIdsToUse(const gmx::MDLogger  &mdlog,
                     ArrayRef<const int>   compatibleGpuIds,
                     const int             numGpuDevicesDetected,
                     const std::string    &userGpuIdChoicesString);

/*! \brief If the user made a GPU ID task assignment, throw if it
 * could not be parsed, or does not suit the available GPU IDs.
 *
 * The error is given with a suitable descriptive message, which will
 * have context if this check is done after the hardware detection
 * results have been reported to the user. However, note that only the
 * GPUs detected on the master rank are reported, because of the
 * existing limitations of that reporting.
 *
 * \param[in]   gpuIdsToUse                  The GPU IDs previously determined as available to use.
 * \param[in]   userGpuTaskAssignmentString  The GPU task assignment from the user.
 *
 * \throws  std::bad_alloc     If out of memory
 *          InvalidInputError  If the assigned GPUs are not valid
 * \returns  The validated user GPU task assignment
 */
std::vector<int>
validateUserGpuTaskAssignment(ArrayRef<const int> gpuIdsToUse,
                              const std::string  &userGpuTaskAssignmentString);

} // namespace

#endif
