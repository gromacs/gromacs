/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_HARDWARE_HWINFO_H
#define GMX_HARDWARE_HWINFO_H

#include <memory>
#include <string>
#include <vector>

#include "gromacs/hardware/device_management.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{
enum class GpuAwareMpiStatus : int;
class CpuInfo;
class HardwareTopology;
} // namespace gmx
struct DeviceInformation;

/* Hardware information structure with CPU and GPU information.
 * It is initialized by gmx_detect_hardware().
 * NOTE: this structure may only contain structures that are
 *       valid over the whole process (i.e. must be able to
 *       be shared among all threads) */
struct gmx_hw_info_t
{
    gmx_hw_info_t(std::unique_ptr<gmx::CpuInfo>          theCpuInfo,
                  std::unique_ptr<gmx::HardwareTopology> theHardwareTopology);
    ~gmx_hw_info_t();

    /* Data for our local physical node */

    std::unique_ptr<gmx::CpuInfo>          cpuInfo; /* Information about CPU capabilities */
    std::unique_ptr<gmx::HardwareTopology> hardwareTopology; /* Information about hardware topology */
    std::vector<std::unique_ptr<DeviceInformation>> deviceInfoList; /* Information about GPUs detected on this physical node */


    /* Data reduced through MPI over all physical nodes */
    int nphysicalnode;        /* Number of physical nodes */
    int ncore_tot;            /* Sum of #cores over all nodes, can be 0 */
    int ncore_min;            /* Min #cores over all nodes */
    int ncore_max;            /* Max #cores over all nodes */
    int nProcessingUnits_tot; /* Sum of # available processing units over all nodes */
    int nProcessingUnits_min; /* Min # available processing units over all nodes */
    int nProcessingUnits_max; /* Max # available processing units over all nodes */
    int maxThreads_tot;       /* Sum of # recommended threads to start over all nodes */
    int maxThreads_min;       /* Min # recommended threads to start over all nodes */
    int maxThreads_max;       /* Max # recommended threads to start over all nodes */
    int ngpu_compatible_tot;  /* Sum of #GPUs over all nodes */
    int ngpu_compatible_min;  /* Min #GPUs over all nodes */
    int ngpu_compatible_max;  /* Max #GPUs over all nodes */

    int simd_suggest_min; /* Highest SIMD instruction set supported by all ranks */
    int simd_suggest_max; /* Highest SIMD instruction set supported by at least one rank */

    gmx_bool bIdenticalGPUs; /* TRUE if all ranks have the same type(s) and order of GPUs */
    bool     haveAmdZen1Cpu; /* TRUE when at least one CPU in any of the nodes is AMD Zen of the first generation */
    gmx::GpuAwareMpiStatus minGpuAwareMpiStatus; /* Lowest support level for GPU-aware MPI (Supported > Forced > Not supported) across all detected devices */

    //! Container of warning strings to log later when that is possible.
    std::vector<std::string> hardwareDetectionWarnings_;
};


/* The options for the thread affinity setting, default: auto */
enum class ThreadAffinity
{
    Select,
    Auto,
    On,
    Off,
    Count
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
    int nthreads_tot = 0;
    //! Number of thread-MPI threads requested.
    int nthreads_tmpi = 0;
    //! Number of OpenMP threads requested.
    int nthreads_omp = 0;
    //! Number of OpenMP threads to use on PME_only ranks.
    int nthreads_omp_pme = 0;
    //! Thread affinity switch, see enum above.
    ThreadAffinity threadAffinity = ThreadAffinity::Select;
    //! Logical core pinning stride.
    int core_pinning_stride = 0;
    //! Logical core pinning offset.
    int core_pinning_offset = 0;
    //! Empty, or a string provided by the user declaring (unique) GPU IDs available for mdrun to use.
    std::string devicesSelectedByUser;
    //! Empty, or a string provided by the user mapping GPU tasks to devices.
    std::string userGpuTaskAssignment;
    //! Tells whether mdrun is free to choose the total number of threads (by choosing the number of OpenMP and/or thread-MPI threads).
    bool totNumThreadsIsAuto;
};

#endif
