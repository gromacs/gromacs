/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016 by the GROMACS development team.
 * Copyright (c) 2017,2018,2019,2020, by the GROMACS development team, led by
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
 * \brief Defines routine for activating potentially deactivated cores
 * so they can be detected.
 *
 * The use of std::thread makes for brittle interaction with std
 * library headers. Its caller also handles GPU detection and
 * allocation of device-specific data structures. This is more
 * manageable when separated into two distinct translation units.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_hardware
 */
#include "gmxpre.h"

#include "prepare_detection.h"

#include "config.h"

#include <cstdio>

#include <chrono>
#include <thread>
#include <vector>

#include "architecture.h"

#ifdef HAVE_UNISTD_H
#    include <unistd.h> // sysconf()
#endif

namespace gmx
{

/*! \brief Utility that does dummy computing for max 2 seconds to spin up cores
 *
 *  This routine will check the number of cores configured and online
 *  (using sysconf), and the spins doing dummy compute operations for up to
 *  2 seconds, or until all cores have come online. This can be used prior to
 *  hardware detection for platforms that take unused processors offline.
 *
 *  This routine will not throw exceptions. In principle it should be
 *  declared noexcept, but at least icc 19.1 and 21-beta08 with the
 *  libstdc++-7.5 has difficulty implementing a std::vector of
 *  std::thread started with this function when declared noexcept. It
 *  is not clear whether the problem is the compiler or the standard
 *  library. Fortunately, this function is not performance sensitive,
 *  and only runs on platforms other than x86 and POWER (ie ARM),
 *  so the possible overhead introduced by omitting noexcept is not
 *  important.
 */
static void spinUpCore()
{
#if defined(HAVE_SYSCONF) && defined(_SC_NPROCESSORS_CONF) && defined(_SC_NPROCESSORS_ONLN)
    float dummy           = 0.1;
    int   countConfigured = sysconf(_SC_NPROCESSORS_CONF);    // noexcept
    auto  start           = std::chrono::steady_clock::now(); // noexcept

    while (sysconf(_SC_NPROCESSORS_ONLN) < countConfigured
           && std::chrono::steady_clock::now() - start < std::chrono::seconds(2))
    {
        for (int i = 1; i < 10000; i++)
        {
            dummy /= i;
        }
    }

    if (dummy < 0)
    {
        printf("This cannot happen, but prevents loop from being optimized away.");
    }
#endif
}

void hardwareTopologyPrepareDetection()
{
#if defined(HAVE_SYSCONF) && defined(_SC_NPROCESSORS_CONF) \
        && (defined(THREAD_PTHREADS) || defined(THREAD_WINDOWS))

    // Modify this conditional when/if x86 or PowerPC starts to sleep some cores
    if (c_architecture != Architecture::X86 && c_architecture != Architecture::PowerPC)
    {
        int                      countConfigured = sysconf(_SC_NPROCESSORS_CONF);
        std::vector<std::thread> workThreads(countConfigured);

        for (auto& t : workThreads)
        {
            t = std::thread(spinUpCore);
        }

        for (auto& t : workThreads)
        {
            t.join();
        }
    }
#endif
}

} // namespace gmx
