/*
 * This file is part of the GROMACS molecular simulation package,
 * version 4.6
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*
 * define HAVE_RDTSCP to use the serializing rdtscp instruction instead of rdtsc.
 * This is only supported on newer Intel/AMD hardware, but provides better accuracy.
 */

/* check for cycle counters on supported platforms */
#if ((defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__)  || defined(__PGIC__)) && (defined(__i386__) || defined(__x86_64__)))
#define TMPI_CYCLE_COUNT
/* x86 or x86-64 with GCC inline assembly */
typedef unsigned long long tmpi_cycles_t;

static __inline__ tmpi_cycles_t tmpi_cycles_read(void)
{
    /* x86 with GCC inline assembly - pentium TSC register */
    tmpi_cycles_t   cycle;
    unsigned       low,high;

#ifdef HAVE_RDTSCP
    __asm__ __volatile__("rdtscp" : "=a" (low), "=d" (high) :: "ecx" );
#else
    __asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high));
#endif

    cycle = ((unsigned long long)low) | (((unsigned long long)high)<<32);

    return cycle;
}
#elif (defined(__INTEL_COMPILER) && defined(__ia64__))
#define TMPI_CYCLE_COUNT
typedef unsigned long tmpi_cycles_t;
static __inline__ tmpi_cycles_t tmpi_cycles_read(void)
{
    /* Intel compiler on ia64 */
    return __getReg(_IA64_REG_AR_ITC);
}
#elif defined(__GNUC__) && defined(__ia64__)
#define TMPI_CYCLE_COUNT
typedef unsigned long tmpi_cycles_t;
static __inline__ tmpi_cycles_t tmpi_cycles_read(void)
{
    /* ia64 with GCC inline assembly */
    tmpi_cycles_t ret;
    __asm__ __volatile__ ("mov %0=ar.itc" : "=r"(ret));
    return ret;
}
#elif defined(_MSC_VER)
#define TMPI_CYCLE_COUNT
typedef __int64 tmpi_cycles_t;
static __inline tmpi_cycles_t tmpi_cycles_read(void)
{
#ifdef HAVE_RDTSCP
    unsigned int ui;
    return __rdtscp(&ui);
#else
	return __rdtsc();
#endif
}
#endif

