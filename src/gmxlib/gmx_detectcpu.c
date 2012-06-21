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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifdef _MSC_VER
/* MSVC definition for __cpuid() */
#include <intrin.h>
#endif


#include "gmx_detectcpu.h"


const char *
gmx_detectcpu_vendorid_string[GMX_DETECTCPU_NVENDORS] =
{
    "Unknown",
    "GenuineIntel",
    "AuthenticAMD"
};

const char *
gmx_detectcpu_feature_string[GMX_DETECTCPU_NFEATURES] =
{
    "CannotDetect",
    "htt",
    "sse2",
    "sse4.1",
    "rdrand",
    "aes",
    "avx",
    "fma",
    "fma4",
    "xop",
    "avx2",
    "rdtscp"
};

const char *
gmx_detectcpu_acceleration_string[GMX_DETECTCPU_NACCELERATIONS] =
{
    "None",
    "SSE2",
    "SSE4.1",
    "AVX_128_FMA",
    "AVX_256"
};





/* What type of acceleration was compiled in, if any?
 * This is set from Cmake. Note that the SSE2 and SSE4_1 macros are set for
 * AVX too, so it is important that they appear last in the list.
 */
#ifdef GMX_X86_AVX_256
static const
gmx_detectcpu_acceleration_t 
compiled_acc = GMX_DETECTCPU_ACCELERATION_X86_AVX_256;
#elif defined GMX_X86_AVX_128_FMA
static const
gmx_detectcpu_acceleration_t 
compiled_acc = GMX_DETECTCPU_ACCELERATION_X86_AVX_128_FMA;
#elif defined GMX_X86_SSE4_1
static const
gmx_detectcpu_acceleration_t 
compiled_acc = GMX_DETECTCPU_ACCELERATION_X86_SSE4_1;
#elif defined GMX_X86_SSE2
static const
gmx_detectcpu_acceleration_t 
compiled_acc = GMX_DETECTCPU_ACCELERATION_X86_SSE2;
#else
static const
gmx_detectcpu_acceleration_t 
compiled_acc = GMX_DETECTCPU_ACCELERATION_NONE;
#endif

/* Execute CPUID on x86 class CPUs. level sets function to exec, and the
 * contents of register output is returned. See Intel/AMD docs for details.
 */
#if defined (__i386__) || defined (__x86_64__) || defined (_M_IX86) || defined (_M_X64)
/* Currently CPUID is only supported (1) if we can use an instruction on MSVC, or (2)
 * if the compiler handles GNU-style inline assembly.
 */
#if (defined GMX_X86_GCC_INLINE_ASM || defined _MSC_VER)
#define GMX_X86_HAVE_CPUID
static int
execute_cpuid_x86(unsigned int level,
                  unsigned int * eax,
                  unsigned int * ebx,
                  unsigned int * ecx,
                  unsigned int * edx)
{
    unsigned int _eax,_ebx,_ecx,_edx;
    int rc;

#ifdef _MSC_VER
    int CPUInfo[4];

    /* MSVC */
    __cpuid(CPUInfo,level);

    _eax=CPUInfo[0];
    _ebx=CPUInfo[1];
    _ecx=CPUInfo[2];
    _edx=CPUInfo[3];

    rc = 0;

#else
    /* for now this means GMX_X86_GCC_INLINE_ASM should be defined,
     * but there might be more options added in the future.
     */
    /* tested on 32 & 64 GCC, and Intel icc. */
#if defined (__x86_64__) || defined (_M_X64)
    __asm__("pushl %%rbx      \n\t"
            "cpuid            \n\t"
            "movl %%ebx, %1   \n\t"
            "popl %%rbx       \n\t"
            : "=a"(_eax), "=r"(_ebx), "=c"(_ecx), "=d"(_edx) : "0"(level));
#else
    __asm__("pushl %%ebx      \n\t"
            "cpuid            \n\t"
            "movl %%ebx, %1   \n\t"
            "popl %%ebx       \n\t"
            : "=a"(_eax), "=r"(_ebx), "=c"(_ecx), "=d"(_edx) : "0"(level));
#endif
    
    rc = 0;
#endif
    /* If you end up having a compiler that really doesn't understand this and
     * you can't fix it, create a separate ifdef and set the results to:
     *
     * _eax=_ebx=_ecx=_edx=0;
     * rc = -1;
     *
     * However, this will lose you ALL Gromacs x86 acceleration, so you want to
     * try really hard before giving up!
     */

    *eax = _eax;
    *ebx = _ebx;
    *ecx = _ecx;
    *edx = _edx;

    return rc;
}
#endif /* GMX_X86_GCC_INLINE_ASM or _MSC_VER */
#endif /* architecture is x86 */


/* Identify CPU features common to Intel & AMD - mainly brand string,
 * version and some features. Vendor has already been detected outside this.
 */
static int
detectcpu_common_x86(gmx_detectcpu_t *              data)
{
    int                       fn,max_stdfn,max_extfn;
    unsigned int              eax,ebx,ecx,edx;
    char                      str[GMX_DETECTCPU_STRLEN];
    char *                    p;

#ifdef GMX_X86_HAVE_CPUID
    /* Find largest standard/extended function input value */
    execute_cpuid_x86(0x0,&eax,&ebx,&ecx,&edx);
    max_stdfn = eax;
    execute_cpuid_x86(0x80000000,&eax,&ebx,&ecx,&edx);
    max_extfn = eax;

    p = str;
    if(max_extfn>=0x80000005)
    {
        /* Get CPU brand string */
        for(fn=0x80000002;fn<0x80000005;fn++)
        {
            execute_cpuid_x86(fn,&eax,&ebx,&ecx,&edx);
            memcpy(p,&eax,4);
            memcpy(p+4,&ebx,4);
            memcpy(p+8,&ecx,4);
            memcpy(p+12,&edx,4);
            p+=16;
        }
        *p='\0';

        /* Remove empty initial space */
        p = str;
        while(isspace(*(p)))
        {
            p++;
        }
    }
    else
    {
        *p='\0';
    }
    strncpy(data->brand,p,GMX_DETECTCPU_STRLEN);

    /* Find basic CPU properties */
    if(max_stdfn>=1)
    {
        execute_cpuid_x86(1,&eax,&ebx,&ecx,&edx);

        data->family   = ((eax & 0x0FF00000) >> 20) + ((eax & 0x00000F00) >> 8);
        /* Note that extended model should be shifted left 4, so only shift right 12 iso 16. */
        data->model    = ((eax & 0x000F0000) >> 12) + ((eax & 0x000000F0) >> 4);
        data->stepping = (eax & 0x0000000F);

        /* Feature flags common to AMD and intel */
        data->feature[GMX_DETECTCPU_FEATURE_X86_FMA]     = (ecx & (1 << 12)) != 0;
        data->feature[GMX_DETECTCPU_FEATURE_X86_SSE4_1]  = (ecx & (1 << 19)) != 0;
        data->feature[GMX_DETECTCPU_FEATURE_X86_AES]     = (ecx & (1 << 25)) != 0;
        data->feature[GMX_DETECTCPU_FEATURE_X86_AVX]     = (ecx & (1 << 28)) != 0;
        data->feature[GMX_DETECTCPU_FEATURE_X86_RDRAND]  = (ecx & (1 << 30)) != 0;

        data->feature[GMX_DETECTCPU_FEATURE_X86_SSE2]    = (edx & (1 << 26)) != 0;
        data->feature[GMX_DETECTCPU_FEATURE_X86_HTT]     = (edx & (1 << 28)) != 0;
    }

    if(max_extfn>=0x80000001)
    {
        execute_cpuid_x86(0x80000001,&eax,&ebx,&ecx,&edx);
        data->feature[GMX_DETECTCPU_FEATURE_X86_RDTSCP]  = (edx & (1 << 27)) != 0;
    }

#else
    /* No CPUID present */
    strncpy(data->brand,"Unknown CPU brand",GMX_DETECTCPU_STRLEN);
    data->family   = 0;
    data->model    = 0;
    data->stepping = 0;
#endif

    return 0;
}

/* Detection of AMD-specific CPU features */
static int
detectcpu_amd(gmx_detectcpu_t *              data)
{
    int                       max_stdfn,max_extfn;
    unsigned int              eax,ebx,ecx,edx;

    detectcpu_common_x86(data);

#ifdef GMX_X86_HAVE_CPUID
    execute_cpuid_x86(0x0,&eax,&ebx,&ecx,&edx);
    max_stdfn = eax;

    execute_cpuid_x86(0x80000000,&eax,&ebx,&ecx,&edx);
    max_extfn = eax;

    if(max_extfn>=0x80000001)
    {
        execute_cpuid_x86(0x80000001,&eax,&ebx,&ecx,&edx);

        data->feature[GMX_DETECTCPU_FEATURE_X86_XOP]     = (ecx & (1 << 11)) != 0;
        data->feature[GMX_DETECTCPU_FEATURE_X86_FMA4]    = (ecx & (1 << 16)) != 0;
    }
#endif

    return 0;
}

/* Detection of Intel-specific CPU features */
static int
detectcpu_intel(gmx_detectcpu_t *              data)
{
    int                       max_stdfn;
    unsigned int              eax,ebx,ecx,edx;

    detectcpu_common_x86(data);

#ifdef GMX_X86_HAVE_CPUID
    execute_cpuid_x86(0x0,&eax,&ebx,&ecx,&edx);
    max_stdfn = eax;

    if(max_stdfn>=7)
    {
        execute_cpuid_x86(0x7,&eax,&ebx,&ecx,&edx);
        data->feature[GMX_DETECTCPU_FEATURE_X86_AVX2]    = (ebx & (1 << 5)) != 0;
    }

#endif

    return 0;
}

/* Try to find the vendor of the current CPU, so we know what specific
 * detection routine to call.
 */
static gmx_detectcpu_vendorid_t
detectcpu_vendor(void)
{
    gmx_detectcpu_vendorid_t   i,vendor;
    /* Register data used on x86 */
    unsigned int               eax,ebx,ecx,edx;
    char                       vendorstring[13];

    /* Set default first */
    vendor = GMX_DETECTCPU_VENDOR_UNKNOWN;

#ifdef GMX_X86_HAVE_CPUID
    execute_cpuid_x86(0,&eax,&ebx,&ecx,&edx);

    memcpy(vendorstring,&ebx,4);
    memcpy(vendorstring+4,&edx,4);
    memcpy(vendorstring+8,&ecx,4);

    vendorstring[12]='\0';

    for(i=GMX_DETECTCPU_VENDOR_UNKNOWN;i<GMX_DETECTCPU_NVENDORS;i++)
    {
        if(!strncmp(vendorstring,gmx_detectcpu_vendorid_string[i],12))
        {
            vendor = i;
        }
    }
#endif

    return vendor;
}

int
gmx_detectcpu                   (gmx_detectcpu_t *              data)
{
    int i;

    for(i=0;i<GMX_DETECTCPU_NFEATURES;i++)
    {
        data->feature[i]=0;
    }
    
    data->vendorid = detectcpu_vendor();

    switch(data->vendorid)
    {
        case GMX_DETECTCPU_VENDOR_INTEL:
            detectcpu_intel(data);
            break;
        case GMX_DETECTCPU_VENDOR_AMD:
            detectcpu_amd(data);
            break;
        default:
            /* Could not find vendor */
            strncpy(data->brand,"Unknown CPU brand",GMX_DETECTCPU_STRLEN);
            data->family         = 0;
            data->model          = 0;
            data->stepping       = 0;

            for(i=0;i<GMX_DETECTCPU_NFEATURES;i++)
            {
                data->feature[i]=0;
            }
            data->feature[GMX_DETECTCPU_FEATURE_CANNOTDETECT] = 1;
            break;
    }

    return 0;
}




int
gmx_detectcpu_formatstring       (gmx_detectcpu_t              data,
                                  char *                        str,
                                  int                           n)
{
    int c;
    int i;

#ifdef _MSC_VER
    _snprintf(str,n,
              "Vendor: %s\n"
              "Brand:  %s\n"
              "Family: %2d  Model: %2d  Stepping: %2d\n"
              "Features:",
              gmx_detectcpu_vendorid_string[data.vendorid],
              data.brand,
              data.family,data.model,data.stepping);
#else
    snprintf(str,n,
             "Vendor: %s\n"
             "Brand:  %s\n"
             "Family: %2d  Model: %2d  Stepping: %2d\n"
             "Features:",
             gmx_detectcpu_vendorid_string[data.vendorid],
             data.brand,
             data.family,data.model,data.stepping);
#endif

    str[n-1] = '\0';
    c = strlen(str);
    n   -= c;
    str += c;

    for(i=0;i<GMX_DETECTCPU_NFEATURES;i++)
    {
        if(data.feature[i]==1)
        {
#ifdef _MSC_VER
            _snprintf(str,n," %s",gmx_detectcpu_feature_string[i]);
#else
            snprintf(str,n," %s",gmx_detectcpu_feature_string[i]);
#endif
            str[n-1] = '\0';
            c = strlen(str);
            n   -= c;
            str += c;
        }
    }
#ifdef _MSC_VER
    _snprintf(str,n,"\n");
#else
    snprintf(str,n,"\n");
#endif
    str[n-1] = '\0';

    return 0;
}



int
gmx_detectcpu_suggest_acceleration  (gmx_detectcpu_t                 data,
                                     gmx_detectcpu_acceleration_t *  acc)
{
    gmx_detectcpu_acceleration_t tmpacc;

    tmpacc = GMX_DETECTCPU_ACCELERATION_NONE;

    if(data.vendorid==GMX_DETECTCPU_VENDOR_INTEL)
    {
        if(data.feature[GMX_DETECTCPU_FEATURE_X86_AVX]==1)
        {
            tmpacc = GMX_DETECTCPU_ACCELERATION_X86_AVX_256;
        }
        else if(data.feature[GMX_DETECTCPU_FEATURE_X86_SSE4_1]==1)
        {
            tmpacc = GMX_DETECTCPU_ACCELERATION_X86_SSE4_1;
        }
        else if(data.feature[GMX_DETECTCPU_FEATURE_X86_SSE2]==1)
        {
            tmpacc = GMX_DETECTCPU_ACCELERATION_X86_SSE2;
        }
    }
    else if(data.vendorid==GMX_DETECTCPU_VENDOR_AMD)
    {
        if(data.feature[GMX_DETECTCPU_FEATURE_X86_AVX]==1)
        {
            tmpacc = GMX_DETECTCPU_ACCELERATION_X86_AVX_128_FMA;
        }
        else if(data.feature[GMX_DETECTCPU_FEATURE_X86_SSE4_1]==1)
        {
            tmpacc = GMX_DETECTCPU_ACCELERATION_X86_SSE4_1;
        }
        else if(data.feature[GMX_DETECTCPU_FEATURE_X86_SSE2]==1)
        {
            tmpacc = GMX_DETECTCPU_ACCELERATION_X86_SSE2;
        }
    }

    *acc = tmpacc;

    return 0;
}



int
gmx_detectcpu_check_acceleration(gmx_detectcpu_t   data,
                                 FILE *           log)
{
    int                           rc;
    char                          str[1024];
    gmx_detectcpu_acceleration_t  acc;

    gmx_detectcpu_suggest_acceleration(data,&acc);
    rc = (acc != compiled_acc);

    gmx_detectcpu_formatstring(data,str,1023);
    str[1023] = '\0';

    if(log!=NULL)
    {
        fprintf(log,
                "Detecting CPU-specific acceleration. Present hardware specification:\n"
                "%s"
                "Acceleration most likely to fit this hardware: %s\n"
                "Acceleration selected at Gromacs compile time: %s\n\n",
                str,
                gmx_detectcpu_acceleration_string[acc],
                gmx_detectcpu_acceleration_string[compiled_acc]);
    }

    if(rc!=0)
    {
        if(log!=NULL)
        {
            fprintf(log,"WARNING! Binary not matching hardware - you are likely losing performance.\n\n");
        }
        printf("\nWARNING! Binary not matching hardware - you are likely losing performance.\n"
               "Acceleration most likely to fit this hardware: %s\n"
               "Acceleration selected at Gromacs compile time: %s\n\n",
               gmx_detectcpu_acceleration_string[acc],
               gmx_detectcpu_acceleration_string[compiled_acc]);
    }

    return rc;
}




#ifdef GMX_DETECTCPU_STANDALONE
/* Stand-alone program to enable queries of CPU features from Cmake.
 * Note that you need to check inline ASM capabilities before compling and set 
 * -DGMX_X86_GCC_INLINE_ASM for the cpuid instruction to work...
 */
int
main(int argc, char **argv)
{
    gmx_detectcpu_t               data;
    gmx_detectcpu_acceleration_t  acc;
    int                           i,cnt;

    if(argc<2)
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
                "-acceleration  Print suggested Gromacs acceleration.\n"
                ,argv[0]);
        exit(0);
    }

    gmx_detectcpu(&data);

    if(!strncmp(argv[1],"-vendor",3))
    {
        printf("%s\n",gmx_detectcpu_vendorid_string[data.vendorid]);
    }
    else if(!strncmp(argv[1],"-brand",3))
    {
        printf("%s\n",data.brand);
    }
    else if(!strncmp(argv[1],"-family",3))
    {
        printf("%d\n",data.family);
    }
    else if(!strncmp(argv[1],"-model",3))
    {
        printf("%d\n",data.model);
    }
    else if(!strncmp(argv[1],"-stepping",3))
    {
        printf("%d\n",data.stepping);
    }
    else if(!strncmp(argv[1],"-features",3))
    {
        cnt = 0;
        for(i=0;i<GMX_DETECTCPU_NFEATURES;i++)
        {
            if(data.feature[i]==1)
            {
                if(cnt++ > 0)
                {
                    printf(" ");
                }
                printf("%s",gmx_detectcpu_feature_string[i]);
            }
        }
        printf("\n");
    }
    else if(!strncmp(argv[1],"-acceleration",3))
    {
        gmx_detectcpu_suggest_acceleration(data,&acc);
        fprintf(stdout,"%s\n",gmx_detectcpu_acceleration_string[acc]);
    }

    return 0;
}

#endif
