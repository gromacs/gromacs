/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include "gmxpre.h"

#include "cyclecounter.h"

#include "config.h"

#include <ctime>

#ifdef HAVE_SYS_TIME_H
#    include <sys/time.h>
#endif
#ifdef _MSC_VER
#    include <windows.h>
#endif

#include "gromacs/utility/basedefinitions.h"

/*! \brief Calculate number of seconds per cycle tick on host
 *
 *  This routine runs a timer loop to calibrate the number of
 *  seconds per the units returned fro gmx_cycles_read().
 *
 *  \param  sampletime Minimum real sample time. It takes some trial-and-error
 *          to find the correct delay loop size, so the total runtime of
 *          this routine is about twice this time.
 *  \return Number of seconds per cycle unit. If it is not possible to
 *          calculate on this system (for whatever reason) the return value
 *          will be -1, so check that it is positive before using it.
 */
double gmx_cycles_calibrate(double sampletime)
{
    /* On ARM and recent-generation x86-64, we can use the more accurate cycle counters
     * that allow better timing for things that depend on it (e.g. load balancing, profiling).
     */
#if ((defined __aarch64__) \
     && (defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__) || defined(__PGIC__)))
    /* 64-bit ARM cycle counters with GCC inline assembly */
    unsigned long cycles;
    __asm__ __volatile__("mrs %0, cntfrq_el0" : "=r"(cycles));
    /* Only first 32 bits are significant */
    cycles &= 0xFFFFFFFF;
    return 1. / cycles;
    GMX_UNUSED_VALUE(sampletime);
#else
#    if ((defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__) || defined(__PGIC__)) \
         && defined(__x86_64__) && !defined(__ILP32__) && !defined(_CRAYC))
    long gmx_unused tmp;
    int             cpuid1;
    int gmx_unused  cpuid2;
    const int       l0  = 0x0;
    const int       l16 = 0x16;
    gmx_cycles_t    cycles;

    /* cpuid clobbers ebx but it must be restored for -fPIC so save
     * then restore ebx */
    __asm__ volatile(
            "xchg %%rbx, %2\n"
            "cpuid\n"
            "xchg %%rbx, %2\n"
            : "=a"(cpuid1), "=d"(cpuid2), "=r"(tmp)
            : "a"(l0)
            : "ecx", "ebx");
    if (cpuid1 >= 0x16)
    {
        /* This CPU is recent enough so the timer frequency can be directly queried */
        __asm__ volatile(
                "xchg %%rbx, %2\n"
                "cpuid\n"
                "xchg %%rbx, %2\n"
                : "=a"(cpuid1), "=d"(cpuid2), "=r"(tmp)
                : "a"(l16)
                : "ecx", "ebx");
        cycles = static_cast<gmx_cycles_t>(cpuid1) * static_cast<gmx_cycles_t>(1000000);
        return 1. / cycles;
    }
#    endif
#    ifdef _MSC_VER

    /* Windows does not have gettimeofday, but it provides a special
     * routine that returns the cycle counter frequency.
     */
    LARGE_INTEGER i;

    QueryPerformanceFrequency(&i);

    return 1.0 / static_cast<double>(i.QuadPart);
    /* end of MS Windows implementation */

#    elif HAVE_GETTIMEOFDAY

    /*  generic implementation with gettimeofday() */
    struct timeval t1, t2;
    gmx_cycles_t   c1, c2;
    double         timediff, cyclediff;
    double         d = 0.1; /* Dummy variable so we don't optimize away delay loop */

    if (!gmx_cycles_have_counter())
    {
        return -1;
    }

#        if (defined(__alpha__) || defined(__alpha))
    /* Alpha cannot count to more than 4e9, but I don't expect
     * that the architecture will go over 2GHz before it dies, so
     * up to 2.0 seconds of sampling should be safe.
     */
    if (sampletime > 2.0)
    {
        sampletime = 2.0;
    }
#        endif

    /* Start a timing loop. We want this to be largely independent
     * of machine speed, so we need to start with a very small number
     * of iterations and repeat it until we reach the requested time.
     *
     * We call gettimeofday an extra time at the start to avoid cache misses.
     */
    gettimeofday(&t1, nullptr);
    gettimeofday(&t1, nullptr);
    c1 = gmx_cycles_read();

    do
    {
        /* just a delay loop. To avoid optimizing it away, we calculate a number
         * that will underflow to zero in most cases. By conditionally adding it
         * to a result at the end it cannot be removed. n=10000 is arbitrary...
         */
        for (int i = 0; i < 10000; i++)
        {
            d = d / (1.0 + static_cast<double>(i));
        }
        /* Read the time again */
        gettimeofday(&t2, nullptr);
        c2       = gmx_cycles_read();
        timediff = static_cast<double>(t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) * 1e-6;
    } while (timediff < sampletime);

    cyclediff = c2 - c1;

    /* Add a very small result so the delay loop cannot be optimized away */
    if (d < 1e-30)
    {
        timediff += d;
    }

    /* Return seconds per cycle */
    return timediff / cyclediff;

#    else
    /* No timing function available */
    return -1;
    GMX_UNUSED_VALUE(sampletime);
#    endif
#endif
}
