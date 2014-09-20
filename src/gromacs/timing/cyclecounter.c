/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2006 David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "cyclecounter.h"

#include "config.h"

#include <time.h>

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef _MSC_VER
#include <windows.h>
#endif

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
double
gmx_cycles_calibrate(double sampletime)
{
#ifdef _MSC_VER

    /* Windows does not have gettimeofday, but it provides a special
     * routine that returns the cycle counter frequency.
     */
    LARGE_INTEGER i;

    QueryPerformanceFrequency(&i);

    return 1.0/((double) i.QuadPart);
    /* end of MS Windows implementation */

#elif (defined HAVE_GETTIMEOFDAY)

    /*  generic implementation with gettimeofday() */
    struct timeval t1, t2;
    gmx_cycles_t   c1, c2;
    double         timediff, cyclediff;
    double         d = 0.1; /* Dummy variable so we don't optimize away delay loop */
    int            i;

    if (!gmx_cycles_have_counter())
    {
        return -1;
    }

#if (defined(__alpha__) || defined(__alpha))
    /* Alpha cannot count to more than 4e9, but I don't expect
     * that the architecture will go over 2GHz before it dies, so
     * up to 2.0 seconds of sampling should be safe.
     */
    if (sampletime > 2.0)
    {
        sampletime = 2.0;
    }
#endif

    /* Start a timing loop. We want this to be largely independent
     * of machine speed, so we need to start with a very small number
     * of iterations and repeat it until we reach the requested time.
     *
     * We call gettimeofday an extra time at the start to avoid cache misses.
     */
    gettimeofday(&t1, NULL);
    gettimeofday(&t1, NULL);
    c1 = gmx_cycles_read();

    do
    {
        /* just a delay loop. To avoid optimizing it away, we calculate a number
         * that will underflow to zero in most cases. By conditionally adding it
         * to a result at the end it cannot be removed. n=10000 is arbitrary...
         */
        for (i = 0; i < 10000; i++)
        {
            d = d/(1.0+(double)i);
        }
        /* Read the time again */
        gettimeofday(&t2, NULL);
        c2       = gmx_cycles_read();
        timediff = (double)(t2.tv_sec-t1.tv_sec)+
            (double)(t2.tv_usec-t1.tv_usec)*1e-6;
    }
    while (timediff < sampletime);

    cyclediff = c2-c1;

    /* Add a very small result so the delay loop cannot be optimized away */
    if (d < 1e-30)
    {
        timediff += d;
    }

    /* Return seconds per cycle */
    return timediff/cyclediff;

#else
    /* No timing function available */
    return -1;
#endif
}
