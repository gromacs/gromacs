/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2006
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifndef _GMX_CYCLECOUNTER_H_
#define _GMX_CYCLECOUNTER_H_

/** @file gmx_cyclecounter.h
 *
 *  @brief High-resolution timestamp or CPU clock cycle counters.
 *
 *  After reading the current value with gmx_cycles_read() you can add or
 *  subtract these numbers as normal integers of type gmx_cycles_t.
 */

#ifdef _MSC_VER
#include <intrin.h>
#endif

#ifdef __cplusplus
extern "C" 
{
#endif
#if 0
} /* fixes auto-indentation problems */
#endif



/* Minor implementation note: 
 *
 * I like to use these counters in other programs too, so to avoid making
 * it dependent on other Gromacs definitions I use the #ifdef's to set
 * architecture-specific inline macros instead of using gmx_inline from
 * gmx_types.h /Erik 2005-12-10
 */

#if ((defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__) || defined(__PGIC__)) && \
     (defined(__i386__) || defined(__x86_64__)))
/* x86 or x86-64 with GCC inline assembly */
typedef unsigned long long 
gmx_cycles_t;

#elif defined(_MSC_VER)
#include <windows.h>
typedef __int64
gmx_cycles_t;

#elif (defined(__hpux) || defined(__HP_cc)) && defined(__ia64)
/* HP compiler on ia64 */
#include <machine/sys/inline.h>
typedef unsigned long 
gmx_cycles_t;

#elif (defined(__INTEL_COMPILER) || defined(__ECC)) && defined(__ia64__) 
/* Intel compiler on ia64 */
#include <ia64intrin.h>
typedef unsigned long 
gmx_cycles_t;

#elif defined(__GNUC__) && defined(__ia64__)
/* ia64 with GCC inline assembly */
typedef unsigned long 
gmx_cycles_t;

#elif ((defined(__hppa__) || defined(__hppa)) && defined (__GNUC__))
/* HP PA-RISC, inline asm with gcc */
typedef unsigned long 
gmx_cycles_t;

#elif ((defined(__hppa__) || defined(__hppa)) && defined (__hpux))
/* HP PA-RISC, instruction when using HP compiler */
#include <machine/inline.h>
typedef unsigned long 
gmx_cycles_t;

#elif defined(__GNUC__) && defined(__s390__) 
/* S390, taken from FFTW who got it from James Treacy */
typedef unsigned long long 
gmx_cycles_t;

#elif defined(__GNUC__) && defined(__alpha__)
/* gcc inline assembly on alpha CPUs */
typedef unsigned long
gmx_cycles_t;

#elif defined(__GNUC__) && defined(__sparc_v9__)
/* gcc inline assembly on sparc v9 */
typedef unsigned long
gmx_cycles_t;

#elif defined(__DECC) && defined(__alpha) 
/* Digital GEM C compiler on alpha */
#include <c_asm.h>
typedef unsigned long
gmx_cycles_t;

#elif (defined(__sgi) && defined(CLOCK_SGI_CYCLE))
/* Irix compilers on SGI hardware. Get nanoseconds from struct timespec */
typedef unsigned long long
gmx_cycles_t;

#elif (defined(__SVR4) && defined (__SUNPRO_CC))
/* Solaris high-resolution timers */
typedef hrtime_t 
gmx_cycles_t;

#elif defined(__xlC__) && defined (_AIX)
/* AIX compilers */
#include <sys/time.h>
#include <sys/systemcfg.h>
typedef unsigned long long
gmx_cycles_t;

#elif ( ( defined(__GNUC__) || defined(__IBM_GCC_ASM) || defined(__IBM_STDCPP_ASM) ) && \
        ( defined(__powerpc__) || defined(__ppc__) ) )
/* PowerPC using gcc inline assembly (also works on xlc>=7.0 with -qasm=gcc) */
typedef unsigned long long 
gmx_cycles_t;

#elif (defined(__MWERKS__) && (defined(MAC) || defined(macintosh)))
/* Metrowerks on macintosh */
typedef unsigned long long 
gmx_cycles_t;

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
typedef long
gmx_cycles_t;

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
#if ((defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__) || defined(__PGIC__)) && \
     (defined(__i386__) || defined(__x86_64__)))
static __inline__ int gmx_cycles_have_counter(void)
{ 
	/* x86 or x86-64 with GCC inline assembly - pentium TSC register */
	return 1;
}
#elif (defined(_MSC_VER))
static __inline int gmx_cycles_have_counter(void)
{ 
	return 1;
}
#elif (defined(__hpux) || defined(__HP_cc)) && defined(__ia64)
static inline int gmx_cycles_have_counter(void)
{ 
	/* HP compiler on ia64, use special instruction to read ITC */
	return 1;
}
#elif (defined(__INTEL_COMPILER) || defined(__ECC)) && defined(__ia64__) 
static __inline__ int gmx_cycles_have_counter(void)
{ 
	/* Intel compiler on ia64, use special instruction to read ITC */
	return 1;
}
#elif defined(__GNUC__) && defined(__ia64__)
static __inline__ int gmx_cycles_have_counter(void)
{ 
	/* AMD64 with GCC inline assembly - TSC register */
	return 1;
}
#elif ((defined(__hppa__) || defined(__hppa)) && defined (__GNUC__))
static __inline__ int gmx_cycles_have_counter(void)
{ 
	/* HP PA-RISC, inline asm with gcc */
	return 1;
}
#elif ((defined(__hppa__) || defined(__hppa)) && defined (__hpux))
static inline int gmx_cycles_have_counter(void)
{ 
	/* HP PA-RISC, instruction when using HP compiler */
	return 1;
}
#elif defined(__GNUC__) && defined(__s390__) 
static __inline__ int gmx_cycles_have_counter(void)
{ 
    /* S390, taken from FFTW who got it from James Treacy */
	return 1;
}
#elif defined(__GNUC__) && defined(__alpha__)
static __inline__ int gmx_cycles_have_counter(void)
{ 
    /* gcc inline assembly on alpha CPUs */
    return 1;
}
#elif defined(__GNUC__) && defined(__sparc_v9__)
static __inline__ int gmx_cycles_have_counter(void)
{ 
    /* gcc inline assembly on sparc v9 */
    return 1;
}
#elif defined(__DECC) && defined(__alpha) 
static __inline int gmx_cycles_have_counter(void)
{ 
    /* Digital GEM C compiler on alpha */
    return 1;
}
#elif (defined(__sgi) && defined(CLOCK_SGI_CYCLE))
static __inline int gmx_cycles_have_counter(void)
{ 
    /* Irix compilers on SGI hardware */
    return 1;
}
#elif (defined(__SVR4) && defined (__SUNPRO_CC))
static inline int gmx_cycles_have_counter(void)
{ 
    /* Solaris high-resolution timers */
    return 1;
}
#elif defined(__xlC__) && defined (_AIX)
static inline int gmx_cycles_have_counter(void)
{ 
    /* AIX compilers */
    return 1;
}
#elif ( ( defined(__GNUC__) || defined(__IBM_GCC_ASM) || defined(__IBM_STDCPP_ASM) ) && \
        ( defined(__powerpc__) || defined(__ppc__) ) )
static __inline__ int gmx_cycles_have_counter(void)
{ 
    /* PowerPC using gcc inline assembly (and xlc>=7.0 with -qasm=gcc) */
    return 1;
}
#elif (defined(__MWERKS__) && (defined(MAC) || defined(macintosh)))
static __inline__ int gmx_cycles_have_counter(void)
{ 
    /* Metrowerks on macintosh */
    return 1;
}
#else
static int gmx_cycles_have_counter(void)
{ 
    /* No cycle counter that we know of on this system */
    return 0;
}
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
#if ((defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__) || defined(__PGIC__)) && \
     (defined(__i386__) || defined(__x86_64__)))
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{ 
    /* x86 with GCC inline assembly - pentium TSC register */
    gmx_cycles_t   cycle;
    unsigned       low,high;
    
    __asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high));
    
    cycle = ((unsigned long long)low) | (((unsigned long long)high)<<32); 
    
    return cycle;
}
#elif defined(_MSC_VER)
static __inline gmx_cycles_t gmx_cycles_read(void)
{ 
	return __rdtsc();
}
#elif (defined(__hpux) || defined(__HP_cc)) && defined(__ia64)
static inline gmx_cycles_t gmx_cycles_read(void)
{ 
    /* HP compiler on ia64 */
    gmx_cycles_t ret;
    ret = _Asm_mov_from_ar (_AREG_ITC);
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
    __asm__ __volatile__ ("mov %0=ar.itc" : "=r"(ret));
    return ret;
}
#elif ((defined(__hppa__) || defined(__hppa)) && defined (__GNUC__))
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{ 
    /* HP PA-RISC, inline asm with gcc */
    gmx_cycles_t ret;
    __asm__ __volatile__("mfctl 16, %0": "=r" (ret));
    /* no input, nothing else clobbered */
    return ret;
}
#elif ((defined(__hppa__) || defined(__hppa)) && defined (__hpux))
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
    __asm__("stck 0(%0)" : : "a" (&(cycle)) : "memory", "cc");
    return cycle;
}
#elif defined(__GNUC__) && defined(__alpha__)
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{ 
    /* gcc inline assembly on alpha CPUs */
    unsigned long cycle;
    __asm__ __volatile__ ("rpcc %0" : "=r"(cycle));
    return (cycle & 0xFFFFFFFF);
}
#elif defined(__GNUC__) && defined(__sparc_v9__)
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{ 
    /* gcc inline assembly on sparc v9 */
    unsigned long ret;
    __asm__("rd %%tick, %0" : "=r" (ret));
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
    return ((unsigned long long)t.tv_sec)*1000000000+
        (unsigned long long)t.tv_nsec;
}
#elif (defined(__SVR4) && defined (__SUNPRO_CC))
static inline gmx_cycles_t gmx_cycles_read(void)
{ 
    /* Solaris high-resolution timers */
    return gethrtime();
}
#elif defined(__xlC__) && defined (_AIX)
static inline gmx_cycles_t gmx_cycles_read(void)
{ 
    /* AIX compilers. Inline the calculation instead of using library functions */
    timebasestruct_t t1;
    read_real_time(&t1, TIMEBASE_SZ);
    /* POWER returns real time (seconds + nanoseconds),
     * POWER_PC returns high/low 32 bits of a counter.
     */
    if(t1.flag==RTC_POWER_PC) 
    {
        return ((gmx_cycles_t)t1.tb_high)<<32 | (gmx_cycles_t)t1.tb_low;
    }
    else
    {
        return ((gmx_cycles_t)t1.tb_high)*1000000000+(gmx_cycles_t)t1.tb_low;
    }
}
#elif ( ( defined(__GNUC__) || defined(__IBM_GCC_ASM) || defined(__IBM_STDCPP_ASM) ) && \
        ( defined(__powerpc__) || defined(__ppc__) ) )
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{ 
    /* PowerPC using gcc inline assembly (and xlC>=7.0 with -qasm=gcc) */
    unsigned long low, high1, high2;
    do
    {
        __asm__ __volatile__ ("mftbu %0" : "=r" (high1) : );
        __asm__ __volatile__ ("mftb %0" : "=r" (low) : );
        __asm__ __volatile__ ("mftbu %0" : "=r" (high2) : );
    } 
    while (high1 != high2);
    
    return (((gmx_cycles_t)high2) << 32) | (gmx_cycles_t)low;
} 
#elif (defined(__MWERKS__) && (defined(MAC) || defined(macintosh)))
static __inline__ gmx_cycles_t gmx_cycles_read(void)
{ 
    /* Metrowerks on macintosh */
    unsigned int long low, high1, high2;
    do
    {
        __asm__ __volatile__ ("mftbu %0" : "=r" (high1) : );
        __asm__ __volatile__ ("mftb %0" : "=r" (low) : );
        __asm__ __volatile__ ("mftbu %0" : "=r" (high2) : );
    } 
    while (high1 != high2);
    
    return (((gmx_cycles_t)high2) << 32) | (gmx_cycles_t)low;  
}
#else
static gmx_cycles_t gmx_cycles_read(void)
{ 
    return 0;
}
#endif








/*! \brief Calculate number of seconds per cycle tick on host
*
*  This routine runs a timer loop to calibrate the number of
*  seconds per the units returned from gmx_cycles_difference()
*
*  To calculate the time used, call gmx_cycles_read() twice,
*  and then use this routine to calculate the difference as a double
*  precision floating-point number.
*
*  \param  sampletime Minimum number of seconds to sample. 
*          One second should give you a reasonably accurate calibration.
*  \return Number of seconds per cycle unit. If it is not possible to
*          calculate on this system (for whatever reason) the return value
*          will be -1, so check that it is positive before using it.
*/
double 
gmx_cycles_calibrate(double sampletime);


#ifdef __cplusplus
}
#endif



#endif /* _GMX_CYCLECOUNTER_H_ */

