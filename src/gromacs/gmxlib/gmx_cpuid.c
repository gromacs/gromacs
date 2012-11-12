/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 * This file is part of GROMACS.
 * Copyright (c) 2012-
 *
 * Written by the Gromacs development team under coordination of
 * David van der Spoel, Berk Hess, and Erik Lindahl.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 *
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_SCHED_H
#define _GNU_SOURCE
#include <sched.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifdef _MSC_VER
/* MSVC definition for __cpuid() */
#include <intrin.h>
#endif
#ifdef HAVE_UNISTD_H
/* sysconf() definition */
#include <unistd.h>
#endif




#include "gmx_cpuid.h"


/* Global constant character strings corresponding to our enumerated types */
const char *
gmx_cpuid_vendor_string[GMX_CPUID_NVENDORS] = {
    "CannotDetect",
    "Unknown",
    "GenuineIntel",
    "AuthenticAMD"
};

const char *
gmx_cpuid_feature_string[GMX_CPUID_NFEATURES] = {
    "CannotDetect",
    "aes",
    "apic",
    "avx",
    "avx2",
    "clfsh",
    "cmov",
    "cx8",
    "cx16",
    "f16c",
    "fma",
    "fma4",
    "htt",
    "lahf_lm",
    "misalignsse",
    "mmx",
    "msr",
    "nonstop_tsc",
    "pcid",
    "pclmuldq",
    "pdcm",
    "pdpe1gb",
    "popcnt",
    "pse",
    "rdrnd",
    "rdtscp",
    "sse2",
    "sse3",
    "sse4a",
    "sse4.1",
    "sse4.2",
    "ssse3",
    "tdt",
    "x2apic",
    "xop"
};

const char *
gmx_cpuid_acceleration_string[GMX_CPUID_NACCELERATIONS] = {
    "CannotDetect",
    "None",
    "SSE2",
    "SSE4.1",
    "AVX_128_FMA",
    "AVX_256"
};

/* Max length of brand string */
#define GMX_CPUID_BRAND_MAXLEN 256


/* Contents of the abstract datatype */
struct gmx_cpuid
{
    enum gmx_cpuid_vendor vendor;
    char brand[GMX_CPUID_BRAND_MAXLEN];
    int family;
    int model;
    int stepping;
    /* Not using gmx_bool here, since this file must be possible to compile without simple.h */
    char feature[GMX_CPUID_NFEATURES];
};


/* Simple routines to access the data structure. The initialization routine is
 * further down since that needs to call other static routines in this file.
 */
enum gmx_cpuid_vendor
gmx_cpuid_vendor            (gmx_cpuid_t cpuid)
{
    return cpuid->vendor;
}


const char *
gmx_cpuid_brand             (gmx_cpuid_t cpuid)
{
    return cpuid->brand;
}

int
gmx_cpuid_family            (gmx_cpuid_t cpuid)
{
    return cpuid->family;
}

int
gmx_cpuid_model             (gmx_cpuid_t cpuid)
{
    return cpuid->model;
}

int
gmx_cpuid_stepping          (gmx_cpuid_t cpuid)
{
    return cpuid->stepping;
}

int
gmx_cpuid_feature           (gmx_cpuid_t            cpuid,
                             enum gmx_cpuid_feature feature)
{
    return (cpuid->feature[feature] != 0);
}




/* What type of acceleration was compiled in, if any?
 * This is set from Cmake. Note that the SSE2 and SSE4_1 macros are set for
 * AVX too, so it is important that they appear last in the list.
 */
#ifdef GMX_X86_AVX_256
static const
enum gmx_cpuid_acceleration
    compiled_acc = GMX_CPUID_ACCELERATION_X86_AVX_256;
#elif defined GMX_X86_AVX_128_FMA
static const
enum gmx_cpuid_acceleration
    compiled_acc = GMX_CPUID_ACCELERATION_X86_AVX_128_FMA;
#elif defined GMX_X86_SSE4_1
static const
enum gmx_cpuid_acceleration
    compiled_acc = GMX_CPUID_ACCELERATION_X86_SSE4_1;
#elif defined GMX_X86_SSE2
static const
enum gmx_cpuid_acceleration
    compiled_acc = GMX_CPUID_ACCELERATION_X86_SSE2;
#else
static const
enum gmx_cpuid_acceleration
    compiled_acc = GMX_CPUID_ACCELERATION_NONE;
#endif


/* Currently CPUID is only supported (1) if we can use an instruction on MSVC, or (2)
 * if the compiler handles GNU-style inline assembly.
 */
#if defined (__i386__) || defined (__x86_64__) || defined (_M_IX86) || defined (_M_X64)

/* Execute CPUID on x86 class CPUs. level sets function to exec, and the
 * contents of register output is returned. See Intel/AMD docs for details.
 *
 * This version supports extended information where we can also have an input
 * value in the ecx register. This is ignored for most levels, but some of them
 * (e.g. level 0xB on Intel) use it.
 */
static int
execute_x86cpuid(unsigned int   level,
                 unsigned int   ecxval,
                 unsigned int * eax,
                 unsigned int * ebx,
                 unsigned int * ecx,
                 unsigned int * edx)
{
    int rc = 0;

#if (defined _MSC_VER)
    int CPUInfo[4];

#if (_MSC_VER > 1500) || (_MSC_VER == 1500 & _MSC_FULL_VER >= 150030729)
    /* MSVC 9.0 SP1 or later */
    __cpuidex(CPUInfo, level, ecxval);
    rc = 0;
#else
    __cpuid(CPUInfo, level);
    /* Set an error code if the user wanted a non-zero ecxval, since we did not have cpuidex */
    rc   = (ecxval > 0) ? -1 : 0;
#endif
    *eax = CPUInfo[0];
    *ebx = CPUInfo[1];
    *ecx = CPUInfo[2];
    *edx = CPUInfo[3];

#elif (defined GMX_X86_GCC_INLINE_ASM)
    /* for now this means GMX_X86_GCC_INLINE_ASM should be defined,
     * but there might be more options added in the future.
     */
    *eax = level;
    *ecx = ecxval;
    *ebx = 0;
    *edx = 0;
#if defined(__i386__) && defined(__PIC__)
    /* Avoid clobbering the global offset table in 32-bit pic code (ebx register) */
    __asm__ __volatile__ ("xchgl %%ebx, %1  \n\t"
                          "cpuid            \n\t"
                          "xchgl %%ebx, %1  \n\t"
                          : "+a" (*eax), "+r" (*ebx), "+c" (*ecx), "+d" (*edx));
#else
    /* i386 without PIC, or x86-64. Things are easy and we can clobber any reg we want :-) */
    __asm__ __volatile__ ("cpuid            \n\t"
                          : "+a" (*eax), "+b" (*ebx), "+c" (*ecx), "+d" (*edx));
#endif
    rc = 0;
#else
    /* Death and horror!
     * Apparently this is an x86 platform where we don't know how to call cpuid.
     *
     * This is REALLY bad, since we will lose all Gromacs acceleration.
     */
    *eax = 0;
    *ebx = 0;
    *ecx = 0;
    *edx = 0;

    rc   = -1;
#endif
    return rc;
}
#endif /* architecture is x86 */


/* Identify CPU features common to Intel & AMD - mainly brand string,
 * version and some features. Vendor has already been detected outside this.
 */
static int
cpuid_check_common_x86(gmx_cpuid_t cpuid)
{
    int          fn, max_stdfn, max_extfn;
    unsigned int eax, ebx, ecx, edx;
    char         str[GMX_CPUID_BRAND_MAXLEN];
    char       * p;

    /* Find largest standard/extended function input value */
    execute_x86cpuid(0x0, 0, &eax, &ebx, &ecx, &edx);
    max_stdfn = eax;
    execute_x86cpuid(0x80000000, 0, &eax, &ebx, &ecx, &edx);
    max_extfn = eax;

    p         = str;
    if(max_extfn >= 0x80000005)
    {
        /* Get CPU brand string */
        for(fn = 0x80000002; fn < 0x80000005; fn++)
        {
            execute_x86cpuid(fn, 0, &eax, &ebx, &ecx, &edx);
            memcpy(p, &eax, 4);
            memcpy(p+4, &ebx, 4);
            memcpy(p+8, &ecx, 4);
            memcpy(p+12, &edx, 4);
            p += 16;
        }
        *p = '\0';

        /* Remove empty initial space */
        p = str;
        while(isspace(*(p)))
        {
            p++;
        }
        strncpy(cpuid->brand, p, GMX_CPUID_BRAND_MAXLEN);
    }
    else
    {
        strncpy(cpuid->brand, "Unknown CPU brand", GMX_CPUID_BRAND_MAXLEN);
    }

    /* Find basic CPU properties */
    if(max_stdfn >= 1)
    {
        execute_x86cpuid(0x1, 0, &eax, &ebx, &ecx, &edx);

        cpuid->family   = ((eax & 0x0FF00000) >> 20) + ((eax & 0x00000F00) >> 8);
        /* Note that extended model should be shifted left 4, so only shift right 12 iso 16. */
        cpuid->model    = ((eax & 0x000F0000) >> 12) + ((eax & 0x000000F0) >> 4);
        cpuid->stepping = (eax & 0x0000000F);

        /* Feature flags common to AMD and intel */
        cpuid->feature[GMX_CPUID_FEATURE_X86_SSE3]     = (ecx & (1 << 0))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_PCLMULDQ] = (ecx & (1 << 1))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_SSSE3]    = (ecx & (1 << 9))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_FMA]      = (ecx & (1 << 12)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_CX16]     = (ecx & (1 << 13)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_SSE4_1]   = (ecx & (1 << 19)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_SSE4_2]   = (ecx & (1 << 20)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_POPCNT]   = (ecx & (1 << 23)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_AES]      = (ecx & (1 << 25)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_AVX]      = (ecx & (1 << 28)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_F16C]     = (ecx & (1 << 29)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_RDRND]    = (ecx & (1 << 30)) != 0;

        cpuid->feature[GMX_CPUID_FEATURE_X86_PSE]      = (edx & (1 << 3))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_MSR]      = (edx & (1 << 5))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_CX8]      = (edx & (1 << 8))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_APIC]     = (edx & (1 << 9))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_CMOV]     = (edx & (1 << 15)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_CLFSH]    = (edx & (1 << 19)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_MMX]      = (edx & (1 << 23)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_SSE2]     = (edx & (1 << 26)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_HTT]      = (edx & (1 << 28)) != 0;
    }
    else
    {
        cpuid->family   = -1;
        cpuid->model    = -1;
        cpuid->stepping = -1;
    }

    if(max_extfn >= 0x80000001)
    {
        execute_x86cpuid(0x80000001, 0, &eax, &ebx, &ecx, &edx);
        cpuid->feature[GMX_CPUID_FEATURE_X86_LAHF_LM] = (ecx & (1 << 0))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_PDPE1GB] = (edx & (1 << 26)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_RDTSCP]  = (edx & (1 << 27)) != 0;
    }

    if(max_extfn >= 0x80000007)
    {
        execute_x86cpuid(0x80000007, 0, &eax, &ebx, &ecx, &edx);
        cpuid->feature[GMX_CPUID_FEATURE_X86_NONSTOP_TSC] = (edx & (1 << 8))  != 0;
    }

    return 0;
}

/* Detection of AMD-specific CPU features */
static int
cpuid_check_amd_x86(gmx_cpuid_t cpuid)
{
    int          max_stdfn, max_extfn;
    unsigned int eax, ebx, ecx, edx;

    cpuid_check_common_x86(cpuid);

    execute_x86cpuid(0x0, 0, &eax, &ebx, &ecx, &edx);
    max_stdfn = eax;

    execute_x86cpuid(0x80000000, 0, &eax, &ebx, &ecx, &edx);
    max_extfn = eax;

    if(max_extfn >= 0x80000001)
    {
        execute_x86cpuid(0x80000001, 0, &eax, &ebx, &ecx, &edx);

        cpuid->feature[GMX_CPUID_FEATURE_X86_SSE4A]       = (ecx & (1 << 6))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_MISALIGNSSE] = (ecx & (1 << 7))  != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_XOP]         = (ecx & (1 << 11)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_FMA4]        = (ecx & (1 << 16)) != 0;
    }

    return 0;
}

/* Detection of Intel-specific CPU features */
static int
cpuid_check_intel_x86(gmx_cpuid_t cpuid)
{
    unsigned int max_stdfn, max_extfn;
    unsigned int eax, ebx, ecx, edx;
    unsigned int i;
    unsigned int max_logical_cores, max_physical_cores;

    cpuid_check_common_x86(cpuid);

    execute_x86cpuid(0x0, 0, &eax, &ebx, &ecx, &edx);
    max_stdfn = eax;

    execute_x86cpuid(0x80000000, 0, &eax, &ebx, &ecx, &edx);
    max_extfn = eax;

    if(max_stdfn >= 1)
    {
        execute_x86cpuid(0x1, 0, &eax, &ebx, &ecx, &edx);
        cpuid->feature[GMX_CPUID_FEATURE_X86_PDCM]   = (ecx & (1 << 15)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_PCID]   = (ecx & (1 << 17)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_X2APIC] = (ecx & (1 << 21)) != 0;
        cpuid->feature[GMX_CPUID_FEATURE_X86_TDT]    = (ecx & (1 << 24)) != 0;
    }

    if(max_stdfn >= 7)
    {
        execute_x86cpuid(0x7, 0, &eax, &ebx, &ecx, &edx);
        cpuid->feature[GMX_CPUID_FEATURE_X86_AVX2] = (ebx & (1 << 5))  != 0;
    }

    /* Check whether Hyper-Threading is enabled, not only supported */
    if(cpuid->feature[GMX_CPUID_FEATURE_X86_HTT] && max_stdfn >= 4)
    {
        execute_x86cpuid(0x1, 0, &eax, &ebx, &ecx, &edx);
        max_logical_cores  = (ebx >> 16) & 0x0FF;
        execute_x86cpuid(0x4, 0, &eax, &ebx, &ecx, &edx);
        max_physical_cores = ((eax >> 26) & 0x3F) + 1;

        /* Clear HTT flag if we only have 1 logical core per physical */
        if(max_logical_cores/max_physical_cores < 2)
        {
            cpuid->feature[GMX_CPUID_FEATURE_X86_HTT] = 0;
        }
    }
    return 0;
}

/* Try to find the vendor of the current CPU, so we know what specific
 * detection routine to call.
 */
static enum gmx_cpuid_vendor
cpuid_check_vendor(void)
{
    enum gmx_cpuid_vendor i, vendor;
    /* Register data used on x86 */
    unsigned int          eax, ebx, ecx, edx;
    char                  vendorstring[13];

    /* Set default first */
    vendor = GMX_CPUID_VENDOR_UNKNOWN;

    execute_x86cpuid(0x0, 0, &eax, &ebx, &ecx, &edx);

    memcpy(vendorstring, &ebx, 4);
    memcpy(vendorstring+4, &edx, 4);
    memcpy(vendorstring+8, &ecx, 4);

    vendorstring[12] = '\0';

    for(i = GMX_CPUID_VENDOR_UNKNOWN; i < GMX_CPUID_NVENDORS; i++)
    {
        if(!strncmp(vendorstring, gmx_cpuid_vendor_string[i], 12))
        {
            vendor = i;
        }
    }

    return vendor;
}




int
gmx_cpuid_init               (gmx_cpuid_t * pcpuid)
{
    gmx_cpuid_t cpuid;
    int         i;

    cpuid   = malloc(sizeof(*cpuid));

    *pcpuid = cpuid;

    for(i = 0; i < GMX_CPUID_NFEATURES; i++)
    {
        cpuid->feature[i] = 0;
    }

    cpuid->vendor = cpuid_check_vendor();

    switch(cpuid->vendor)
    {
    case GMX_CPUID_VENDOR_INTEL:
        cpuid_check_intel_x86(cpuid);
        break;
    case GMX_CPUID_VENDOR_AMD:
        cpuid_check_amd_x86(cpuid);
        break;
    default:
        /* Could not find vendor */
        strncpy(cpuid->brand, "Unknown CPU brand", GMX_CPUID_BRAND_MAXLEN);
        cpuid->family   = 0;
        cpuid->model    = 0;
        cpuid->stepping = 0;

        for(i = 0; i < GMX_CPUID_NFEATURES; i++)
        {
            cpuid->feature[i] = 0;
        }
        cpuid->feature[GMX_CPUID_FEATURE_CANNOTDETECT] = 1;
        break;
    }

    return 0;
}



void
gmx_cpuid_done               (gmx_cpuid_t cpuid)
{
    free(cpuid);
}


int
gmx_cpuid_formatstring       (gmx_cpuid_t cpuid,
                              char      * str,
                              int         n)
{
    int                    c;
    int                    i;
    enum gmx_cpuid_feature feature;

#ifdef _MSC_VER
    _snprintf(str, n,
              "Vendor: %s\n"
              "Brand:  %s\n"
              "Family: %2d  Model: %2d  Stepping: %2d\n"
              "Features:",
              gmx_cpuid_vendor_string[gmx_cpuid_vendor(cpuid)],
              gmx_cpuid_brand(cpuid),
              gmx_cpuid_family(cpuid), gmx_cpuid_model(cpuid), gmx_cpuid_stepping(cpuid));
#else
    snprintf(str, n,
             "Vendor: %s\n"
             "Brand:  %s\n"
             "Family: %2d  Model: %2d  Stepping: %2d\n"
             "Features:",
             gmx_cpuid_vendor_string[gmx_cpuid_vendor(cpuid)],
             gmx_cpuid_brand(cpuid),
             gmx_cpuid_family(cpuid), gmx_cpuid_model(cpuid), gmx_cpuid_stepping(cpuid));
#endif

    str[n-1] = '\0';
    c        = strlen(str);
    n       -= c;
    str     += c;

    for(feature = GMX_CPUID_FEATURE_CANNOTDETECT; feature < GMX_CPUID_NFEATURES; feature++)
    {
        if(gmx_cpuid_feature(cpuid, feature) == 1)
        {
#ifdef _MSC_VER
            _snprintf(str, n, " %s", gmx_cpuid_feature_string[feature]);
#else
            snprintf(str, n, " %s", gmx_cpuid_feature_string[feature]);
#endif
            str[n-1] = '\0';
            c        = strlen(str);
            n       -= c;
            str     += c;
        }
    }
#ifdef _MSC_VER
    _snprintf(str, n, "\n");
#else
    snprintf(str, n, "\n");
#endif
    str[n-1] = '\0';

    return 0;
}



enum gmx_cpuid_acceleration
gmx_cpuid_acceleration_suggest  (gmx_cpuid_t cpuid)
{
    enum gmx_cpuid_acceleration tmpacc;

    tmpacc = GMX_CPUID_ACCELERATION_NONE;

    if(gmx_cpuid_vendor(cpuid) == GMX_CPUID_VENDOR_INTEL)
    {
        if(gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_AVX))
        {
            tmpacc = GMX_CPUID_ACCELERATION_X86_AVX_256;
        }
        else if(gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_SSE4_1))
        {
            tmpacc = GMX_CPUID_ACCELERATION_X86_SSE4_1;
        }
        else if(gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_SSE2))
        {
            tmpacc = GMX_CPUID_ACCELERATION_X86_SSE2;
        }
    }
    else if(gmx_cpuid_vendor(cpuid) == GMX_CPUID_VENDOR_AMD)
    {
        if(gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_AVX))
        {
            tmpacc = GMX_CPUID_ACCELERATION_X86_AVX_128_FMA;
        }
        else if(gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_SSE4_1))
        {
            tmpacc = GMX_CPUID_ACCELERATION_X86_SSE4_1;
        }
        else if(gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_SSE2))
        {
            tmpacc = GMX_CPUID_ACCELERATION_X86_SSE2;
        }
    }

    return tmpacc;
}



int
gmx_cpuid_acceleration_check(gmx_cpuid_t cpuid,
                             FILE      * log)
{
    int  rc;
    char str[1024];
    enum gmx_cpuid_acceleration acc;

    acc = gmx_cpuid_acceleration_suggest(cpuid);

    rc  = (acc != compiled_acc);

    gmx_cpuid_formatstring(cpuid, str, 1023);
    str[1023] = '\0';

    if(log != NULL)
    {
        fprintf(log,
                "\nDetecting CPU-specific acceleration.\nPresent hardware specification:\n"
                "%s"
                "Acceleration most likely to fit this hardware: %s\n"
                "Acceleration selected at GROMACS compile time: %s\n\n",
                str,
                gmx_cpuid_acceleration_string[acc],
                gmx_cpuid_acceleration_string[compiled_acc]);
    }

    if(rc != 0)
    {
        if(log != NULL)
        {
            fprintf(log, "WARNING! Binary not matching hardware - you are likely losing performance.\n\n");
        }
        printf("\nWARNING! Binary not matching hardware - you are likely losing performance.\n"
               "Acceleration most likely to fit this hardware: %s\n"
               "Acceleration selected at GROMACS compile time: %s\n\n",
               gmx_cpuid_acceleration_string[acc],
               gmx_cpuid_acceleration_string[compiled_acc]);
    }

    return rc;
}


enum gmx_cpuid_x86_smt
gmx_cpuid_x86_smt(gmx_cpuid_t cpuid)
{

#if (defined HAVE_SCHED_H && defined HAVE_SCHED_SETAFFINITY && defined HAVE_SYSCONF && defined __linux__)
    int          i;
    int          nproc;
    cpu_set_t    cpuset, save_cpuset;
    int        * apic_id;
    unsigned int eax, ebx, ecx, edx;
    int          core_shift_bits;
    int          smt_found;

    if( gmx_cpuid_vendor(cpuid) != GMX_CPUID_VENDOR_INTEL ||
        gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_HTT) == 0)
    {
        return GMX_CPUID_X86_SMT_DISABLED;
    }

    /* Check cpuid max standard function */
    execute_x86cpuid(0x0, 0, &eax, &ebx, &ecx, &edx);

    /* Early CPUs that do not support function 11 do not support SMT either */
    if(eax < 0xB)
    {
        return GMX_CPUID_X86_SMT_DISABLED;
    }

    /* If we got here, it is a modern Intel CPU that supports detection, as does our OS */

    /* How many processors? */
    nproc   = sysconf(_SC_NPROCESSORS_ONLN);

    apic_id = malloc(sizeof(int)*nproc);

    sched_getaffinity(0, sizeof(cpu_set_t), &save_cpuset);

    /* Get x2APIC ID from each hardware thread */
    CPU_ZERO(&cpuset);
    for(i = 0; i < nproc; i++)
    {
        CPU_SET(i, &cpuset);
        sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
        execute_x86cpuid(0xB, 0, &eax, &ebx, &ecx, &edx);
        apic_id[i] = edx;
        CPU_CLR(i, &cpuset);
    }
    /* Reset affinity to the value it had when calling this routine */
    sched_setaffinity(0, sizeof(cpu_set_t), &save_cpuset);

    core_shift_bits = eax & 0x1F;

    /* Check if there is any other APIC id that is identical to [0], apart from
     * the hardware thread bit.
     */
    smt_found = 0;
    for(i = 1; i < nproc && smt_found == 0; i++)
    {
        smt_found = (apic_id[i]>>core_shift_bits == apic_id[0] >> core_shift_bits);
    }

    free(apic_id);

    if(smt_found == 1)
    {
        return GMX_CPUID_X86_SMT_ENABLED;
    }
    else
    {
        return GMX_CPUID_X86_SMT_DISABLED;
    }
#else
    /* Do the trivial stuff first. If Hyper-Threading isn't even supported it
     * cannot be enabled, no matter what OS detection we use!
     */
    if(0 == gmx_cpuid_feature(cpuid, GMX_CPUID_FEATURE_X86_HTT))
    {
        return GMX_CPUID_X86_SMT_DISABLED;
    }
    else
    {
        return GMX_CPUID_X86_SMT_CANNOTDETECT;
    }
#endif
}




#ifdef GMX_CPUID_STANDALONE
/* Stand-alone program to enable queries of CPU features from Cmake.
 * Note that you need to check inline ASM capabilities before compiling and set
 * -DGMX_X86_GCC_INLINE_ASM for the cpuid instruction to work...
 */
int
main(int argc, char **argv)
{
    gmx_cpuid_t                 cpuid;
    enum gmx_cpuid_acceleration acc;
    int i, cnt;

    if(argc < 2)
    {
        fprintf(stdout,
                "Usage:\n\n%s [flags]\n\n"
                "Available flags:\n"
                "-vendor        Print CPU vendor.\n"
                "-brand         Print CPU brand string.\n"
                "-family        Print CPU family version.\n"
                "-model         Print CPU model version.\n"
                "-stepping      Print CPU stepping version.\n"
                "-features      Print CPU feature flags.\n"
                "-acceleration  Print suggested GROMACS acceleration.\n",
                argv[0]);
        exit(0);
    }

    gmx_cpuid_init(&cpuid);

    if(!strncmp(argv[1], "-vendor", 3))
    {
        printf("%s\n", gmx_cpuid_vendor_string[cpuid->vendor]);
    }
    else if(!strncmp(argv[1], "-brand", 3))
    {
        printf("%s\n", cpuid->brand);
    }
    else if(!strncmp(argv[1], "-family", 3))
    {
        printf("%d\n", cpuid->family);
    }
    else if(!strncmp(argv[1], "-model", 3))
    {
        printf("%d\n", cpuid->model);
    }
    else if(!strncmp(argv[1], "-stepping", 3))
    {
        printf("%d\n", cpuid->stepping);
    }
    else if(!strncmp(argv[1], "-features", 3))
    {
        cnt = 0;
        for(i = 0; i < GMX_CPUID_NFEATURES; i++)
        {
            if(cpuid->feature[i] == 1)
            {
                if(cnt++ > 0)
                {
                    printf(" ");
                }
                printf("%s", gmx_cpuid_feature_string[i]);
            }
        }
        printf("\n");
    }
    else if(!strncmp(argv[1], "-acceleration", 3))
    {
        acc = gmx_cpuid_acceleration_suggest(cpuid);
        fprintf(stdout, "%s\n", gmx_cpuid_acceleration_string[acc]);
    }

    gmx_cpuid_done(cpuid);


    return 0;
}

#endif
