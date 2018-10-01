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
#ifndef GMX_HARDWARE_HWINFO_H
#define GMX_HARDWARE_HWINFO_H

#include <string>
#include <vector>

#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{
class CpuInfo;
class HardwareTopology;
} // namespace

/* Hardware information structure with CPU and GPU information.
 * It is initialized by gmx_detect_hardware().
 * NOTE: this structure may only contain structures that are globally valid
 *       (i.e. must be able to be shared among all threads) */
struct gmx_hw_info_t
{
    /* Data for our local physical node */
    struct gmx_gpu_info_t gpu_info;                /* Information about GPUs detected in the system */

    int                   nthreads_hw_avail;       /* Number of hardware threads available; this number
                                                      is based on the number of CPUs reported as available
                                                      by the OS at the time of detection. */

    const gmx::CpuInfo *         cpuInfo;          /* Information about CPU capabilities */
    const gmx::HardwareTopology *hardwareTopology; /* Information about hardware topology */

    /* Data reduced through MPI over all physical nodes */
    int                 nphysicalnode;       /* Number of physical nodes */
    int                 ncore_tot;           /* Sum of #cores over all nodes, can be 0 */
    int                 ncore_min;           /* Min #cores over all nodes */
    int                 ncore_max;           /* Max #cores over all nodes */
    int                 nhwthread_tot;       /* Sum of #hwthreads over all nodes */
    int                 nhwthread_min;       /* Min #hwthreads over all nodes */
    int                 nhwthread_max;       /* Max #hwthreads over all nodes */
    int                 ngpu_compatible_tot; /* Sum of #GPUs over all nodes */
    int                 ngpu_compatible_min; /* Min #GPUs over all nodes */
    int                 ngpu_compatible_max; /* Max #GPUs over all nodes */

    int                 simd_suggest_min;    /* Highest SIMD instruction set supported by all ranks */
    int                 simd_suggest_max;    /* Highest SIMD instruction set supported by at least one rank */

    gmx_bool            bIdenticalGPUs;      /* TRUE if all ranks have the same type(s) and order of GPUs */
    bool                haveAmdZenCpu;       /* TRUE when at least one CPU in any of the nodes is AMD Zen */
};


/* The options for the thread affinity setting, default: auto */
enum {
    threadaffSEL, threadaffAUTO, threadaffON, threadaffOFF, threadaffNR
};

/*! \internal \brief Threading and GPU options, can be set automatically or by the user
 *
 * \todo During mdrunner(), if the user has left any of these values
 * at their defaults (which tends to mean "choose automatically"),
 * then those values are over-written with the result of such
 * automation. This creates problems for the subsequent code in
 * knowing what was done, why, and reporting correctly to the
 * user. Find a way to improve this.
 */
struct gmx_hw_opt_t
{
    //! Total number of threads requested (thread-MPI + OpenMP).
    int           nthreads_tot = 0;
    //! Number of thread-MPI threads requested.
    int           nthreads_tmpi = 0;
    //! Number of OpenMP threads requested.
    int           nthreads_omp = 0;
    //! Number of OpenMP threads to use on PME_only ranks.
    int           nthreads_omp_pme = 0;
    //! Thread affinity switch, see enum above.
    int           thread_affinity = threadaffSEL;
    //! Logical core pinning stride.
    int           core_pinning_stride = 0;
    //! Logical core pinning offset.
    int           core_pinning_offset = 0;
    //! Empty, or a string provided by the user declaring (unique) GPU IDs available for mdrun to use.
    std::string   gpuIdsAvailable = "";
    //! Empty, or a string provided by the user mapping GPU tasks to devices.
    std::string   userGpuTaskAssignment = "";
    //! Tells whether mdrun is free to choose the total number of threads (by choosing the number of OpenMP and/or thread-MPI threads).
    bool          totNumThreadsIsAuto;
};

#endif
