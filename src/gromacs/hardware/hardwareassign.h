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

#include <string>
#include <vector>

#include "gromacs/utility/basedefinitions.h"

struct gmx_gpu_info_t;
struct gmx_hw_opt_t;
struct t_commrec;

namespace gmx
{

/*! \brief Parse a GPU assignment string into digits
 *
 * \param[in]   gpuTaskAssignment  String like "013" or "0,1,3" typically
 *                                 supplied by the user to mdrun -gpu_id.
 *
 * \returns  A vector of integer GPU ids, like {0, 1, 3}.
 *
 * \throws   std::bad_alloc     If out of memory.
 *           InvalidInputError  If an invalid character is found (ie not a digit or ',').
 */
std::vector<int> parseGpuTaskAssignment(const std::string &gpuTaskAssignment);

/*! \brief Assign PP ranks to valid GPU IDs.
 *
 * Will return a validated mapping from PP ranks (ie tasks that can
 * run on GPUs) to the device IDs of compatible GPUs on their node.
 * This will be from any non-empty assignment in hw_opt, otherwise a
 * default automated mapping is generated.
 *
 * Note that PME-only ranks have always ignored mdrun -gpu_id, so do
 * not attempt to validate -gpu_id. They should continue this behaviour
 * until PME tasks can use GPUs.
 *
 * \param[in]     rankCanUseGpu          Whether this rank can execute a task on a GPU.
 * \param[in]     cr                     Communication record.
 * \param[in]     gpu_info               Information detected about GPUs, including compatibility.
 * \param[in]     hw_opt                 Parallelisation options, including any user-specified GPU task assignment.
 *
 * \returns  A valid GPU selection.
 */
std::vector<int> mapPpRanksToGpus(bool                    rankCanUseGpu,
                                  const t_commrec        *cr,
                                  const gmx_gpu_info_t   &gpu_info,
                                  const gmx_hw_opt_t     &hw_opt);

} // namespace

#endif
