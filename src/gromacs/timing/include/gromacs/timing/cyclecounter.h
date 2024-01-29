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
/*! \libinternal \file
 * \brief
 * High-resolution timestamp or CPU clock cycle counters.
 *
 * After reading the current value with gmx_cycles_read() you can add or
 * subtract these numbers as normal integers of type gmx_cycles_t.
 *
 * \inlibraryapi
 */
#ifndef GMX_TIMING_CYCLECOUNTER_H
#define GMX_TIMING_CYCLECOUNTER_H

/*
 * Define GMX_USE_RDTSCP=1 to use the serializing rdtscp instruction instead of rdtsc.
 * This is supported on essentially all Intel/AMD hardware still in use, and provides better accuracy.
 */
#include "config.h"

#ifdef _MSC_VER
#    include <intrin.h>
#endif

#if ((defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__) || defined(__PGIC__)) \
     && (defined(__i386__) || defined(__x86_64__)))
/* x86 or x86-64 with GCC inline assembly */
typedef unsigned long long gmx_cycles_t;

#elif ((defined __aarch64__) \
       && (defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__) || defined(__PGIC__)))
/* 64-bit ARM cycle counters with GCC inline assembly */
typedef unsigned long long     gmx_cycles_t;

#elif defined(__ARM_ARCH_7A__) && defined(__GNUC__)
/* Armv7A can provide 64-bit cycles by returning two registers */
typedef unsigned long long     gmx_cycles_t;

#elif defined(_MSC_VER)
#    include <windows.h>
typedef __int64              gmx_cycles_t;

#elif (defined(__hpux) || defined(__HP_cc)) && defined(__ia64)
/* HP compiler on ia64 */
#    include <machine/sys/inline.h>
typedef unsigned long      gmx_cycles_t;

#elif (defined(__INTEL_COMPILER) || defined(__ECC)) && defined(__ia64__)
/* Intel compiler on ia64 */
#    include <ia64intrin.h>
typedef unsigned long          gmx_cycles_t;

#elif defined(__GNUC__) && defined(__ia64__)
/* ia64 with GCC inline assembly */
typedef unsigned long          gmx_cycles_t;

#elif ((defined(__hppa__) || defined(__hppa)) && defined(__GNUC__))
/* HP PA-RISC, inline asm with gcc */
typedef unsigned long          gmx_cycles_t;

#elif ((defined(__hppa__) || defined(__hppa)) && defined(__hpux))
/* HP PA-RISC, instruction when using HP compiler */
#    include <machine/inline.h>
typedef unsigned long      gmx_cycles_t;

#elif defined(__GNUC__) && defined(__s390__)
/* S390, taken from FFTW who got it from James Treacy */
typedef unsigned long long     gmx_cycles_t;

#elif defined(__GNUC__) && defined(__alpha__)
/* gcc inline assembly on alpha CPUs */
typedef unsigned long          gmx_cycles_t;

#elif defined(__GNUC__) && defined(__sparc_v9__)
/* gcc inline assembly on sparc v9 */
typedef unsigned long          gmx_cycles_t;

#elif defined(__DECC) && defined(__alpha)
/* Digital GEM C compiler on alpha */
#    include <c_asm.h>
typedef unsigned long        gmx_cycles_t;

#elif (defined(__sgi) && defined(CLOCK_SGI_CYCLE))
/* Irix compilers on SGI hardware. Get nanoseconds from struct timespec */
typedef unsigned long long   gmx_cycles_t;

#elif (defined(__SVR4) && defined(__SUNPRO_CC))
/* Solaris high-resolution timers */
typedef hrtime_t           gmx_cycles_t;

#elif defined(__xlC__) && defined(_AIX)
/* AIX compilers */
#    include <sys/systemcfg.h>
#    include <sys/time.h>
typedef unsigned long long gmx_cycles_t;

#elif ((defined(__GNUC__) || defined(__IBM_GCC_ASM) || defined(__IBM_STDCPP_ASM)) \
       && (defined(__powerpc__) || defined(__ppc__)))
/* PowerPC using gcc inline assembly (also works on xlc>=7.0 with -qasm=gcc) */
typedef unsigned long long     gmx_cycles_t;

#elif (defined(__MWERKS__) && (defined(MAC) || defined(macintosh)))
/* Metrowerks on macintosh */
typedef unsigned long long     gmx_cycles_t;

#elif defined(__sun) && defined(__sparcv9)

typedef unsigned long gmx_cycles_t;

#elif defined(__loongarch__) && defined(__GNUC__)
typedef unsigned long long gmx_cycles_t;

#else
/*! \brief Integer-like datatype for cycle counter values
 *
 *  Depending on your system this will usually be something like long long,
 *  or a special cycle datatype from the system header files. It is NOT
 *  necessarily real processor cycles - many systems count in nanoseconds
 *  or a special external time register at fixed frequency (not the CPU freq.)
 *
 *  You can subtract or add gmx_cycle_t types just as normal integers, and if
 *  you run the calibration routine you can also multiply it with a factor to
 *  translate the cycle data to seconds.
 */
typedef long                 gmx_cycles_t;

#endif

/*! \brief Read CPU cycle counter
 *
 *  This routine returns an abstract datatype containing a
 *  cycle counter timestamp.
 *
 *  \return Opaque data corresponding to a cycle reading.
 *
 *  Please note that on most systems it takes several cycles
 *  to read and return the cycle counters. If you are measuring
 *  small intervals, you can compensate for this time by calling
 *  the routine twice and calculating what the difference is.
 *  Subtract this from your other measurements to get an accurate result.
 *
 *  Use gmx_cycles_difference() to get a real number corresponding to
 *  the difference between two gmx_cycles_t values returned from this
 *  routine.
 */
#if ((defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__) || defined(__PGIC__)) \
     && (defined(__i386__) || defined(__x86_64__)) && !defined(_CRAYC))
static __inline__ gmx_cycles_t gmx_cycles_read()
{
    /* x86 with GCC inline assembly - pentium TSC register */
    unsigned low, high;

#    if GMX_USE_RDTSCP
    __asm__ __volatile__("rdtscp" : "=a"(low), "=d"(high)::"ecx");
#    else
    __asm__ __volatile__("rdtsc" : "=a"(low), "=d"(high));
#    endif
    const gmx_cycles_t c_low  = low;
    const gmx_cycles_t c_high = high;
    return c_low | c_high << 32;
}
#elif ((defined __aarch64__) \
       && (defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__) || defined(__PGIC__)))
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{
    /* 64-bit ARM cycle counters with GCC inline assembly */
    gmx_cycles_t cycle;
    __asm__ __volatile__("mrs %0, cntvct_el0" : "=r"(cycle));

    return cycle;
}
#elif defined(__ARM_ARCH_7A__) && defined(__GNUC__)
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{
    unsigned int cycles_lo, cycles_hi;
    asm volatile("mrrc p15, 1, %0, %1, c14" : "=r"(cycles_lo), "=r"(cycles_hi));
    return ((gmx_cycles_t)cycles_lo) | (((gmx_cycles_t)cycles_hi) << 32);
}
#elif defined(_MSC_VER)
static __inline gmx_cycles_t gmx_cycles_read(void)
{
#    ifdef _M_ARM
    /* Windows on 64-bit ARM */
    return __rdpmccntr64();
#    else
    /* x86 */
#        if GMX_USE_RDTSCP
    unsigned int ui;
    return __rdtscp(&ui);
#        else
    return __rdtsc();
#        endif
#    endif
}
#elif (defined(__hpux) || defined(__HP_cc)) && defined(__ia64)
static inline gmx_cycles_t gmx_cycles_read(void)
{
    /* HP compiler on ia64 */
    gmx_cycles_t ret;
    ret = _Asm_mov_from_ar(_AREG_ITC);
    return ret;
}
#elif (defined(__INTEL_COMPILER) && defined(__ia64__))
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{
    /* Intel compiler on ia64 */
    return __getReg(_IA64_REG_AR_ITC);
}
#elif defined(__GNUC__) && defined(__ia64__)
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{
    /* ia64 with GCC inline assembly */
    gmx_cycles_t ret;
    __asm__ __volatile__("mov %0=ar.itc" : "=r"(ret));
    return ret;
}
#elif ((defined(__hppa__) || defined(__hppa)) && defined(__GNUC__))
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{
    /* HP PA-RISC, inline asm with gcc */
    gmx_cycles_t ret;
    __asm__ __volatile__("mfctl 16, %0" : "=r"(ret));
    /* no input, nothing else clobbered */
    return ret;
}
#elif ((defined(__hppa__) || defined(__hppa)) && defined(__hpux))
static inline gmx_cycles_t gmx_cycles_read(void)
{
    /* HP PA-RISC, instruction when using HP compiler */
    gmx_cycles_t ret;
    _MFCTL(16, ret);
    return ret;
}
#elif defined(__GNUC__) && defined(__s390__)
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{
    /* S390, taken from FFTW who got it from James Treacy */
    gmx_cycles_t cycle;
    __asm__("stck 0(%0)" : : "a"(&(cycle)) : "memory", "cc");
    return cycle;
}
#elif defined(__GNUC__) && defined(__alpha__)
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{
    /* gcc inline assembly on alpha CPUs */
    unsigned long cycle;
    __asm__ __volatile__("rpcc %0" : "=r"(cycle));
    return (cycle & 0xFFFFFFFF);
}
#elif defined(__GNUC__) && defined(__sparc_v9__)
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{
    /* gcc inline assembly on sparc v9 */
    unsigned long ret;
    __asm__("rd %%tick, %0" : "=r"(ret));
    return ret;
}
#elif defined(__DECC) && defined(__alpha)
static __inline gmx_cycles_t gmx_cycles_read(void)
{
    /* Digital GEM C compiler on alpha */
    unsigned long cycle;
    cycle = asm("rpcc %v0");
    return (cycle & 0xFFFFFFFF);
}
#elif (defined(__sgi) && defined(CLOCK_SGI_CYCLE))
static __inline gmx_cycles_t gmx_cycles_read(void)
{
    /* Irix compilers on SGI hardware */
    struct timespec t;
    clock_gettime(CLOCK_SGI_CYCLE, &t);
    /* Return the number of nanoseconds, so we can subtract/add */
    return ((unsigned long long)t.tv_sec) * 1000000000 + (unsigned long long)t.tv_nsec;
}
#elif (defined(__SVR4) && defined(__SUNPRO_CC))
static inline gmx_cycles_t gmx_cycles_read(void)
{
    /* Solaris high-resolution timers */
    return gethrtime();
}
#elif defined(__xlC__) && defined(_AIX)
static inline gmx_cycles_t gmx_cycles_read(void)
{
    /* AIX compilers. Inline the calculation instead of using library functions */
    timebasestruct_t t1;
    read_real_time(&t1, TIMEBASE_SZ);
    /* POWER returns real time (seconds + nanoseconds),
     * POWER_PC returns high/low 32 bits of a counter.
     */
    if (t1.flag == RTC_POWER_PC)
    {
        return ((gmx_cycles_t)t1.tb_high) << 32 | (gmx_cycles_t)t1.tb_low;
    }
    else
    {
        return ((gmx_cycles_t)t1.tb_high) * 1000000000 + (gmx_cycles_t)t1.tb_low;
    }
}
#elif ((defined(__GNUC__) || defined(__IBM_GCC_ASM) || defined(__IBM_STDCPP_ASM)) \
       && (defined(__powerpc__) || defined(__ppc__)))
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{
    /* PowerPC using gcc inline assembly (and xlC>=7.0 with -qasm=gcc, and clang) */
    unsigned long low, high1, high2;
    do
    {
        // clang 3.7 incorrectly warns that mftb* are
        // deprecated. That's not correct - see
        // https://llvm.org/bugs/show_bug.cgi?id=23680.
        __asm__ __volatile__("mftbu %0" : "=r"(high1) :);
        __asm__ __volatile__("mftb %0" : "=r"(low) :);
        __asm__ __volatile__("mftbu %0" : "=r"(high2) :);
    } while (high1 != high2);

    return (((gmx_cycles_t)high2) << 32) | (gmx_cycles_t)low;
}
#elif (defined(__MWERKS__) && (defined(MAC) || defined(macintosh)))
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{
    /* Metrowerks on macintosh */
    unsigned int long low, high1, high2;
    do
    {
        __asm__ __volatile__("mftbu %0" : "=r"(high1) :);
        __asm__ __volatile__("mftb %0" : "=r"(low) :);
        __asm__ __volatile__("mftbu %0" : "=r"(high2) :);
    } while (high1 != high2);

    return (((gmx_cycles_t)high2) << 32) | (gmx_cycles_t)low;
}
#elif defined(__sun) && defined(__sparcv9)

static __inline__ gmx_cycles_t gmx_cycles_read(void)
{
    gmx_cycles_t ret;
    __asm__ __volatile__("rd %%tick, %0" : "=r"(ret));
    return ret;
}

#elif defined(_CRAYC)
#    include <intrinsics.h>

static __inline gmx_cycles_t gmx_cycles_read(void)
{
    return _rtc();
}
#elif defined __riscv && defined __riscv_xlen && (__riscv_xlen == 32)
static __inline gmx_cycles_t gmx_cycles_read(void)
{
    unsigned int long low, high1, high2;
    do
    {
        asm volatile("csrrs %0, 0xc81, x0" : "=r"(high1));
        asm volatile("csrrs %0, 0xc01, x0" : "=r"(low));
        asm volatile("csrrs %0, 0xc81, x0" : "=r"(high2));
    } while (high1 != high2);

    return (((gmx_cycles_t)high2) << 32) | (gmx_cycles_t)low;
}

#elif defined __riscv && defined __riscv_xlen && (__riscv_xlen == 64)
static __inline gmx_cycles_t gmx_cycles_read(void)
{
    gmx_cycles_t ret;
    asm volatile("csrrs %0, 0xc01, x0" : "=r"(ret));
    return ret;
}

#elif defined __loongarch__ && __loongarch64
static __inline gmx_cycles_t gmx_cycles_read(void)
{
    gmx_cycles_t ret;
    asm volatile("rdtime.d %0, $r0" : "=r"(ret));
    return ret;
}

#else

static gmx_cycles_t gmx_cycles_read(void)
{
    return 0;
}
#endif


/*! \brief Check if high-resolution cycle counters are available
 *
 *  Not all architectures provide any way to read timestep counters
 *  in the CPU, and on some it is broken. Although we refer to it
 *  as cycle counters, it is not necessarily given in units of
 *  cycles.
 *
 *  If you notice that system is missing, implement support for it,
 *  find out how to detect the system during preprocessing, and send us a
 *  patch.
 *
 *  \return 1 if cycle counters are available, 0 if not.
 *
 * \note This functions not need to be in the header for performance
 *       reasons, but it is very important that we get exactly the
 *       same detection as for gmx_cycles_read() routines. If you
 *       compile the library with one compiler, and then use a different
 *       one when later linking to the library it might happen that the
 *       library supports cyclecounters but not the headers, or vice versa.
 */
#if ((defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__) \
      || defined(__PGIC__) || defined(_CRAYC))                                 \
     && (defined(__i386__) || defined(__x86_64__)))
static __inline__ bool gmx_cycles_have_counter()
{
    /* x86 or x86-64 with GCC inline assembly - pentium TSC register */
    return true;
}
#elif ((defined __aarch64__) \
       && (defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__) || defined(__PGIC__)))
static __inline bool gmx_cycles_have_counter(void)
{
    /* 64-bit ARM cycle counters with GCC inline assembly */
    return 1;
}
#elif defined(__ARM_ARCH_7A__) && defined(__GNUC__)
static __inline bool gmx_cycles_have_counter(void)
{
    /* Armv7A can provide 64-bit cycles by returning two registers. However, it will not work unless
     * the performance registers have been made available from user space by a kernel module -
     * otherwise it returns 0.
     */
    gmx_cycles_t c0, c1;

    c0 = gmx_cycles_read();
    c1 = gmx_cycles_read();

    /* if both counters return 0, support is not present */
    return (c0 != 0 || c1 != 0);
}
#elif (defined(_MSC_VER))
static __inline bool gmx_cycles_have_counter(void)
{
    return 1;
}
#elif (defined(__hpux) || defined(__HP_cc)) && defined(__ia64)
static inline bool gmx_cycles_have_counter(void)
{
    /* HP compiler on ia64, use special instruction to read ITC */
    return 1;
}
#elif (defined(__INTEL_COMPILER) || defined(__ECC)) && defined(__ia64__)
static __inline__ bool gmx_cycles_have_counter(void)
{
    /* Intel compiler on ia64, use special instruction to read ITC */
    return 1;
}
#elif defined(__GNUC__) && defined(__ia64__)
static __inline__ bool gmx_cycles_have_counter(void)
{
    /* AMD64 with GCC inline assembly - TSC register */
    return 1;
}
#elif ((defined(__hppa__) || defined(__hppa)) && defined(__GNUC__))
static __inline__ bool gmx_cycles_have_counter(void)
{
    /* HP PA-RISC, inline asm with gcc */
    return 1;
}
#elif ((defined(__hppa__) || defined(__hppa)) && defined(__hpux))
static inline bool gmx_cycles_have_counter(void)
{
    /* HP PA-RISC, instruction when using HP compiler */
    return 1;
}
#elif defined(__GNUC__) && defined(__s390__)
static __inline__ bool gmx_cycles_have_counter(void)
{
    /* S390, taken from FFTW who got it from James Treacy */
    return 1;
}
#elif defined(__GNUC__) && defined(__alpha__)
static __inline__ bool gmx_cycles_have_counter(void)
{
    /* gcc inline assembly on alpha CPUs */
    return 1;
}
#elif defined(__GNUC__) && defined(__sparc_v9__)
static __inline__ bool gmx_cycles_have_counter(void)
{
    /* gcc inline assembly on sparc v9 */
    return 1;
}
#elif defined(__DECC) && defined(__alpha)
static __inline bool gmx_cycles_have_counter(void)
{
    /* Digital GEM C compiler on alpha */
    return 1;
}
#elif (defined(__sgi) && defined(CLOCK_SGI_CYCLE))
static __inline bool gmx_cycles_have_counter(void)
{
    /* Irix compilers on SGI hardware */
    return 1;
}
#elif (defined(__SVR4) && defined(__SUNPRO_CC))
static inline bool gmx_cycles_have_counter(void)
{
    /* Solaris high-resolution timers */
    return 1;
}
#elif defined(__xlC__) && defined(_AIX)
static inline bool gmx_cycles_have_counter(void)
{
    /* AIX compilers */
    return 1;
}
#elif ((defined(__GNUC__) || defined(__IBM_GCC_ASM) || defined(__IBM_STDCPP_ASM)) \
       && (defined(__powerpc__) || defined(__ppc__)))
static __inline__ bool gmx_cycles_have_counter(void)
{
    /* PowerPC using gcc inline assembly (and xlc>=7.0 with -qasm=gcc) */
    return 1;
}
#elif (defined(__MWERKS__) && (defined(MAC) || defined(macintosh)))
static __inline__ bool gmx_cycles_have_counter(void)
{
    /* Metrowerks on macintosh */
    return 1;
}
#elif defined(__sun) && defined(__sparcv9)

static __inline__ bool gmx_cycles_have_counter(void)
{
    /* Solaris on SPARC*/
    return 1;
}
#elif defined __riscv && defined __riscv_xlen && (__riscv_xlen == 32)

static __inline__ bool gmx_cycles_have_counter(void)
{
    /* 32-bit RISC-V */
    return true;
}
#elif defined __riscv && defined __riscv_xlen && (__riscv_xlen == 64)

static __inline__ bool gmx_cycles_have_counter(void)
{
    /* 64-bit RISC-V */
    return true;
}

#elif defined __loongarch__ && __loongarch64

static __inline__ bool gmx_cycles_have_counter(void)
{
    return true;
}

#else
static bool gmx_cycles_have_counter(void)
{
    /* No cycle counter that we know of on this system */
    return 0;
}
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
double gmx_cycles_calibrate(double sampletime);

#endif
