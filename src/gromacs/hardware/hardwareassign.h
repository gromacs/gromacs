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

#include <vector>

#include "gromacs/utility/basedefinitions.h"

struct gmx_gpu_info_t;
struct gmx_gpu_opt_t;
struct t_commrec;

/*! \brief Select the compatible GPUs
 *
 * This function filters gpu_info.gpu_dev for compatible GPUs based
 * on the previously run compatibility tests.
 *
 * \param[in]     gpu_info    Information detected about GPUs, including compatibility
 * \return                    vector of IDs of GPUs already recorded as compatible */
std::vector<int> getCompatibleGpus(const gmx_gpu_info_t &gpu_info);

/*! \brief Map PP ranks to GPU IDs.
 *
 * After this call, gpu_opt->dev_use will contain a validated mapping
 * from PP ranks (ie tasks that can run on GPUs) to the device IDs of
 * compatible GPUs on their node.
 *
 * Note that PME-only ranks have always ignored mdrun -gpu_id, so do
 * not attempt to validate -gpu_id. They should continue this behaviour
 * until PME tasks can use GPUs.
 *
 * \param[in]     rankCanUseGpu  Whether this rank can execute a task on a GPU.
 * \param[in]     cr             Communication record.
 * \param[in]     gpu_info       Information detected about GPUs, including compatibility.
 * \param[in]     userSetGpuIds  Whether the user set the GPU IDs to use in the mapping.
 * \param[inout]  gpu_opt        Holds the mapping to validate, or to fill.
 */
void mapPpRanksToGpus(bool                  rankCanUseGpu,
                      const t_commrec      *cr,
                      const gmx_gpu_info_t &gpu_info,
                      bool                  userSetGpuIds,
                      gmx_gpu_opt_t        *gpu_opt);

#endif
