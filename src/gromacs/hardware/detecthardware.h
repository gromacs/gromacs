/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#ifndef GMX_HARDWARE_DETECTHARDWARE_H
#define GMX_HARDWARE_DETECTHARDWARE_H

#include <cstdio>

#include <string>

#include "gromacs/utility/basedefinitions.h"
// TODO
#include "gromacs/hardware/hardwareassign.h"

struct gmx_gpu_info_t;
struct gmx_gpu_opt_t;
struct gmx_hw_info_t;
struct gmx_hw_opt_t;
struct t_commrec;

namespace gmx
{
class HardwareTopology;
class MDLogger;
}

/*! \brief Run detection, consistency checks, and make available on all ranks.
 *
 * This routine constructs the global hwinfo structure and returns a pointer to
 * it. It will run a preamble before executing cpu and hardware checks, and
 * then run consistency checks afterwards. The results will also be made
 * available on all nodes.
 * Caller is responsible for freeing this pointer.
 */
gmx_hw_info_t *gmx_detect_hardware(const gmx::MDLogger &mdlog,
                                   const t_commrec *cr, gmx_bool bDetectGPUs);

/* Print information about the detected hardware to fplog (if != NULL)
 * and to stderr the master rank.
 */
void gmx_print_detected_hardware(FILE *fplog, const t_commrec *cr,
                                 const gmx::MDLogger &mdlog,
                                 const gmx_hw_info_t *hwinfo);

/*! \brief Helper function for reporting GPU usage information
 * in the mdrun log file
 *
 * \param[in] gpu_info       Information detected about GPUs
 * \param[in] deviceIds      IDs of GPUs assigned to ranks on this node.
 * \param[in] userSetGpuIds  Whether the user selected the GPU ids
 * \param[in] numPpRanks     Number of PP ranks per node
 * \param[in] bPrintHostName Print the hostname in the usage information
 * \return                   String to write to the log file
 * \throws                   std::bad_alloc if out of memory */
std::string
makeGpuUsageReport(const gmx_gpu_info_t &gpu_info,
                   const std::vector<int> &deviceIds,
                   bool                  userSetGpuIds,
                   size_t                numPpRanks,
                   bool                  bPrintHostName);

void gmx_hardware_info_free(gmx_hw_info_t *hwinfo);

/* Return whether the user selected GPU ids */
bool hasUserSetGpuIds(const gmx_gpu_opt_t *gpu_opt);

//! Return whether compatible GPUs were found.
bool compatibleGpusFound(const gmx_gpu_info_t &gpu_info);

/* Parse the GPU ids the user may have passed. */
void gmx_parse_gpu_ids(gmx_gpu_opt_t *gpu_opt);

#endif
