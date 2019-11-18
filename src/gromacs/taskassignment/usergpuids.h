/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019, by the GROMACS development team, led by
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

/*! \brief Parse a GPU ID specifier string into a container describing device IDs exposed to the run.
 *
 * \param[in]   gpuIdString  String like "013" or "0,1,3" typically
 *                           supplied by the user to mdrun -gpu_id.
 *                           Must contain only unique decimal digits, or only decimal
 *                           digits separated by comma delimiters. A terminal
 *                           comma is accceptable (and required to specify a
 *                           single ID that is larger than 9).
 *
 * \returns  A vector of unique GPU IDs.
 *
 * \throws   std::bad_alloc     If out of memory.
 *           InvalidInputError  If an invalid character is found (ie not a digit or ',') or if
 *                              identifiers are duplicated in the specifier list.
 */
std::vector<int> parseUserGpuIdString(const std::string& gpuIdString);

/*! \brief Implement GPU ID selection by returning the available GPU
 * IDs on this physical node that are compatible.
 *
 * If the string supplied by the user is empty, then return the IDs of
 * all compatible GPUs on this physical node. Otherwise, check the
 * user specified compatible GPUs and return their IDs.
 *
 * \param[in]  gpuInfo                Information detected about GPUs on this physical node
 * \param[in]  gpuIdsAvailableString  String like "013" or "0,1,3" typically
 *                                    supplied by the user to mdrun -gpu_id.
 *                                    Must contain only unique decimal digits, or only decimal
 *                                    digits separated by comma delimiters. A terminal
 *                                    comma is accceptable (and required to specify a
 *                                    single ID that is larger than 9).
 *
 * \returns  A vector of unique compatible GPU IDs on this physical node.
 *
 * \throws   std::bad_alloc     If out of memory.
 *           InvalidInputError  If an invalid character is found (ie not a digit or ',') or if
 *                              identifiers are duplicated in the specifier list.
 *           InvalidInputError  If gpuIdsAvailableString specifies GPU IDs that are
 *                              not compatible.
 */
std::vector<int> makeGpuIdsToUse(const gmx_gpu_info_t& gpuInfo, const std::string& gpuIdsAvailableString);

/*! \brief Parse a GPU ID specifier string into a container describing device ID to task mapping.
 *
 * \param[in]   gpuIdString  String like "0011" or "0,0,1,1" typically
 *                           supplied by the user to mdrun -gputasks.
 *                           Must contain only decimal digits, or only decimal
 *                           digits separated by comma delimiters. A terminal
 *                           comma is accceptable (and required to specify a
 *                           single ID that is larger than 9).
 *
 * \returns  A vector of GPU IDs.
 *
 * \throws   std::bad_alloc     If out of memory.
 *           InvalidInputError  If an invalid character is found (ie not a digit or ',').
 */
std::vector<int> parseUserTaskAssignmentString(const std::string& gpuIdString);


/*! \brief Make a vector containing \c numGpuTasks IDs of the IDs found in \c compatibleGpus.
 *
 * \throws  std::bad_alloc          If out of memory
 *
 * \returns A sorted vector of IDs of compatible vectors, whose
 * length matches that of the number of GPU tasks required.
 */
std::vector<int> makeGpuIds(ArrayRef<const int> compatibleGpus, size_t numGpuTasks);

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
std::string makeGpuIdString(const std::vector<int>& gpuIds, int totalNumberOfTasks);

/*! \brief Check that all user-selected GPUs are compatible.
 *
 * Given the \c gpuIds and \c hardwareInfo, throw if
 * any selected GPUs is not compatible.
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
 * \param[in]   gpu_info        Information detected about GPUs
 * \param[in]   compatibleGpus  Vector of GPUs that are compatible
 * \param[in]   gpuIds          The GPU IDs selected by the user.
 *
 * \throws  std::bad_alloc          If out of memory
 *          InconsistentInputError  If the assigned GPUs are not valid
 */
void checkUserGpuIds(const gmx_gpu_info_t&   gpu_info,
                     const std::vector<int>& compatibleGpus,
                     const std::vector<int>& gpuIds);

} // namespace gmx

#endif
