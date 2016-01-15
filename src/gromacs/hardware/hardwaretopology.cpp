/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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

/*! \internal \file
 * \brief
 * Implements gmx::HardwareTopology.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_hardware
 */

#include "gmxpre.h"

#include "hardwaretopology.h"

#include "config.h"

#include <algorithm>
#include <vector>

#include <thread>

#include "gromacs/hardware/cpuinfo.h"

#ifdef HAVE_UNISTD_H
#    include <unistd.h>       // sysconf()
#endif
#if GMX_NATIVE_WINDOWS
#    include <windows.h>      // GetSystemInfo()
#endif

namespace gmx
{

namespace
{


/*! \brief Initialize machine data from basic information in cpuinfo
 *
 *  \param  cpuInfo  CpuInfo object
 *  \param  machine  Machine tree structure where information will be assigned
 *                   if the cpuinfo object contains topology information.
 */
void
parseFromCpuInfo(const gmx::CpuInfo           &cpuInfo,
                 HardwareTopology::Machine *   machine)
{
    if (cpuInfo.logicalProcessors().size() > 0)
    {
        int nSockets   = 0;
        int nCores     = 0;
        int nHwThreads = 0;

        // Copy the logical processor information from cpuinfo
        for (auto &l : cpuInfo.logicalProcessors())
        {
            machine->logicalProcessors.push_back( { l.socket, l.core, l.hwThread } );
            nSockets   = std::max(nSockets, l.socket);
            nCores     = std::max(nCores, l.core);
            nHwThreads = std::max(nHwThreads, l.hwThread);
        }

        // Resize all arrays for sockets/cores/hwthreads properly
        machine->sockets.resize(nSockets + 1);
        for (auto &s : machine->sockets)
        {
            s.cores.resize(nCores + 1);
            for (auto &c : s.cores)
            {
                c.hwThreads.resize(nHwThreads + 1);
            }
        }

        // Fill the logical processor id in the right place
        for (std::size_t i = 0; i < machine->logicalProcessors.size(); i++)
        {
            const HardwareTopology::LogicalProcessor &l = machine->logicalProcessors[i];
            machine->sockets[l.socket].cores[l.core].hwThreads[l.hwThread].logicalProcessorId = static_cast<int>(i);
        }
        machine->logicalProcessorCount = machine->logicalProcessors.size();
    }

}

/*! \brief Try to detect the number of logical processors.
 *
 *  \return The number of hardware processing units, or 0 if it fails.
 */
int
detectLogicalProcessorCount()
{
    // Try to use std::thread::hardware_concurrency() first. This result is only
    // a hint, and it might be 0 if the information is not available.
    // On Apple this will not compile with gcc-4.6, and since it just returns 0 on other
    // platforms too we skip it entirely for gcc < 4.7
#if defined __GNUC__ && (__GNUC__ == 4 && __GNUC_MINOR__ < 7)
    int count = 0;
#else
    int count = std::thread::hardware_concurrency();
#endif

    if (count == 0)
    {
#if GMX_NATIVE_WINDOWS
        // Windows
        SYSTEM_INFO sysinfo;
        GetSystemInfo( &sysinfo );
        count = sysinfo.dwNumberOfProcessors;
#elif defined HAVE_SYSCONF
        // We are probably on Unix. Check if we have the argument to use before executing the call
#    if defined(_SC_NPROCESSORS_CONF)
        count = sysconf(_SC_NPROCESSORS_CONF);
#    elif defined(_SC_NPROC_CONF)
        count = sysconf(_SC_NPROC_CONF);
#    elif defined(_SC_NPROCESSORS_ONLN)
        count = sysconf(_SC_NPROCESSORS_ONLN);
#    elif defined(_SC_NPROC_ONLN)
        count = sysconf(_SC_NPROC_ONLN);
#    endif      // End of check for sysconf argument values

#else
        count = 0; // Neither windows nor Unix, and std::thread_hardware_concurrency() failed.
#endif
    }
    return count;
}

}   // namespace anonymous

// static
HardwareTopology HardwareTopology::detect()
{
    HardwareTopology result;

    CpuInfo          cpuInfo(CpuInfo::detect());

    if (cpuInfo.logicalProcessors().size() > 0)
    {
        // There is topology information in cpuInfo
        parseFromCpuInfo(cpuInfo, &result.machine_);
        result.supportLevel_ = SupportLevel::Basic;
    }
    else
    {
        // No topology information; try to detect the number of logical processors at least
        result.machine_.logicalProcessorCount = detectLogicalProcessorCount();
        if (result.machine_.logicalProcessorCount > 0)
        {
            result.supportLevel_ = SupportLevel::LogicalProcessorCount;
        }
    }

    return result;
}


HardwareTopology::HardwareTopology()
    : supportLevel_(SupportLevel::None)
{
}

} // namespace gmx
