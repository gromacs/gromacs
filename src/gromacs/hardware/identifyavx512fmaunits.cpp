/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief Implements a routine to check the number of AVX512 fma units
 *
 * Just as the CpuInfo code, we need to be able to compile this file in stand-alone mode
 * to set the SIMD acceleration and similar things during CMake configuration.
 */

#ifndef GMX_IDENTIFY_AVX512_FMA_UNITS_STANDALONE
#include "gmxpre.h"
#endif

#include "identifyavx512fmaunits.h"

#ifdef GMX_IDENTIFY_AVX512_FMA_UNITS_STANDALONE
#define AVX_512_FMA_UNIT_DETECTION_COMPILED 1
#else
#include "config.h"
#endif

#include <cstdint>
#include <cstdio>

#include <algorithm>
#include <mutex>

#ifndef GMX_IDENTIFY_AVX512_FMA_UNITS_STANDALONE
#include "gromacs/hardware/cpuinfo.h"
#endif

#include "gromacs/hardware/identifyavx512fmaunits_asm.h"

namespace gmx
{

namespace
{

#if AVX_512_FMA_UNIT_DETECTION_COMPILED
// Use a local routine to read the timestep counter just on x86 to avoid dependence
// on the Gromacs cycle counter module.
uint64_t rdtscp(void)
{
#ifdef MSC_VER
    unsigned int ui;
    return static_cast<uint64_t>(__rdtscp(&ui));
#else
    uint32_t low;
    uint32_t high;

    __asm__ __volatile__("rdtscp" : "=a" (low), "=d" (high) :: "ecx" );
    return (static_cast<uint64_t>(high) << 32) | low;
#endif
}

int
checkDualAvx512FmaUnits()
{
    uint64_t timeFmaAndShuf = 1e9;             // Large value
    uint64_t timeFmaOnly    = 1e9;             // Large value

    // Make sure the CPU is in AVX512 mode by executing a fairly long loop
    executeFmaOnlyLoop(100000);

    // Execute the mixed FMA/shuffle loop three times
    for (int i = 0; i < 3; i++)
    {
        uint64_t start = rdtscp();
        executeFmaAndShuffleLoop(1000);
        uint64_t res = rdtscp() - start;
        timeFmaAndShuf = std::min(timeFmaAndShuf, res);
    }

    // Execute the FMA-only loop three times
    for (int i = 0; i < 3; i++)
    {
        uint64_t start = rdtscp();
        executeFmaOnlyLoop(1000);
        uint64_t res = rdtscp() - start;
        timeFmaOnly = std::min(timeFmaOnly, res);
    }

    // Dummy can never be negative, but by using it in the
    // conditional it cannot be optimized away.
    return (timeFmaAndShuf > 1.5 * timeFmaOnly);
}
#endif  // AVX_512_FMA_UNIT_DETECTION_COMPILED

/*! \brief Mutex to guard the execution of the timing test
 *
 * We only execute the test once, and return the saved result
 * on subsequent calls.
 */
std::mutex initMutex;

}      // namespace anonymous

int
identifyAvx512FmaUnits()
{
    static bool initialized = false;
    static int  result      = false;

    if (!initialized)
    {
        std::lock_guard<std::mutex>  lock(initMutex);

        if (!initialized)
        {
            // For the standalone test binary we assume it will
            // only be executed on AVX512 hardware, but for the
            // library version we check the hardware support.
#ifdef GMX_IDENTIFY_AVX512_FMA_UNITS_STANDALONE
            bool haveAvx512Hardware = true;
#else
            bool haveAvx512Hardware = CpuInfo::detect().feature(CpuInfo::Feature::X86_Avx512F);
#endif

            if (haveAvx512Hardware)
            {
#if AVX_512_FMA_UNIT_DETECTION_COMPILED
                result = checkDualAvx512FmaUnits() ? 2 : 1;
#else
                result = -1; // Cannot run the tests
#endif
            }
            else
            {
                result = 0; // Not AVX-512 hardware
            }
            initialized = true;
        }
    }
    return result;
}

} // namespace gmx

#ifdef GMX_IDENTIFY_AVX512_FMA_UNITS_STANDALONE
int
main()
{
    printf("%d\n", gmx::identifyAvx512FmaUnits());
    return 0;
}
#endif
