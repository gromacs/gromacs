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
#ifndef _GMX_detectcpu_H_
#define _GMX_detectcpu_H_

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif


/* Currently identifiable CPU Vendors */
typedef enum
{
    GMX_DETECTCPU_VENDOR_UNKNOWN = 0,
    GMX_DETECTCPU_VENDOR_INTEL,
    GMX_DETECTCPU_VENDOR_AMD,
    GMX_DETECTCPU_NVENDORS
}
gmx_detectcpu_vendorid_t;

/* Text strings corresponding to CPU vendors */
extern const char *
gmx_detectcpu_vendorid_string[GMX_DETECTCPU_NVENDORS];




/* CPU feature/property list, to be used as indices into the feature array of the
 * gmxDetectCpu_t data structure.
 *
 * Always add entries to the end of this list, just before the last NFEATURES line.
 * To keep the length of this list reasonable, we only add flags referring to
 * features that we actually might have to check/use in Gromacs - feel free to add more.
 *
 * AMD and Intel are unfortunately gradually diverging, so while we can use the
 * same type of intrinsic instruction functions in the source, the resulting binary
 * is frequently not compatible starting from AVX.
 */
typedef enum
{
    GMX_DETECTCPU_FEATURE_CANNOTDETECT = 0,  /* Flag set if we could not detect on this CPU  */
    GMX_DETECTCPU_FEATURE_X86_HTT,           /* Hyperthreading technology                    */
    GMX_DETECTCPU_FEATURE_X86_SSE2,          /* SSE 2                                        */
    GMX_DETECTCPU_FEATURE_X86_SSE4_1,        /* SSE 4.1                                      */
    GMX_DETECTCPU_FEATURE_X86_RDRAND,        /* RDRAND high-quality hardware random numbers  */
    GMX_DETECTCPU_FEATURE_X86_AES,           /* x86 advanced encryption standard accel.      */
    GMX_DETECTCPU_FEATURE_X86_AVX,           /* Advanced vector extensions                   */
    GMX_DETECTCPU_FEATURE_X86_FMA,           /* Fused-multiply add support (mainly for AVX)  */
    GMX_DETECTCPU_FEATURE_X86_FMA4,          /* 4-operand FMA, only on AMD for now           */
    GMX_DETECTCPU_FEATURE_X86_XOP,           /* AMD extended instructions, only AMD for now  */
    GMX_DETECTCPU_FEATURE_X86_AVX2,          /* AVX2 including gather support (not used yet) */
    GMX_DETECTCPU_FEATURE_X86_RDTSCP,        /* Serializing rdtscp instruction available     */
    GMX_DETECTCPU_NFEATURES
}
gmx_detectcpu_feature_t;

/* Text strings for CPU feature indices */
extern const char *
gmx_detectcpu_feature_string[GMX_DETECTCPU_NFEATURES];


/* Currently supported acceleration instruction sets, intrinsics or other similar combinations
 * in Gromacs. There is not always a 1-to-1 correspondence with feature flags; on some AMD
 * hardware we prefer to use 128bit AVX instructions (although 256-bit ones could be executed),
 * and we still havent written the AVX2 kernels.
 */
typedef enum
{
    GMX_DETECTCPU_ACCELERATION_NONE = 0,
    GMX_DETECTCPU_ACCELERATION_X86_SSE2,
    GMX_DETECTCPU_ACCELERATION_X86_SSE4_1,
    GMX_DETECTCPU_ACCELERATION_X86_AVX_128_FMA,
    GMX_DETECTCPU_ACCELERATION_X86_AVX_256,
    GMX_DETECTCPU_NACCELERATIONS
}
gmx_detectcpu_acceleration_t;

/* Text strings for Gromacs acceleration/instruction sets */
extern const char *
gmx_detectcpu_acceleration_string[GMX_DETECTCPU_NACCELERATIONS];



#define GMX_DETECTCPU_STRLEN  64

/* Data structure with CPU detection information. Set by gmxDetectCpu().
 * This is listed in the header for now, since we might want to access it in
 * performance-sensitive part of the code where we don't want function calls.
 */
typedef struct
{
    gmx_detectcpu_vendorid_t  vendorid;
    char                       brand[GMX_DETECTCPU_STRLEN];
    int                        family;
    int                        model;
    int                        stepping;

    char                       feature[GMX_DETECTCPU_NFEATURES];
}
gmx_detectcpu_t;



/* Fill the data structure by using CPU detection instructions.
 * Return 0 on success, 1 if something bad happened.
 */
int
gmx_detectcpu                   (gmx_detectcpu_t *              data);


/* Formats a text string (up to n characters) from the data structure.
 * The output will have max 80 chars between newline characters.
 */
int
gmx_detectcpu_formatstring       (gmx_detectcpu_t                data,
                                  char *                         s,
                                  int                            n);


/* Suggests a suitable gromacs acceleration based on the support in the
 * hardware.
 */
int
gmx_detectcpu_suggest_acceleration  (gmx_detectcpu_t                 data,
                                     gmx_detectcpu_acceleration_t *  acc);

/* Check if this binary was compiled with the same acceleration as we
 * would suggest for the current hardware. Always print stats to the log file
 * if it is non-NULL, and print a warning in stdout if we don't have a match.
 */
int
gmx_detectcpu_check_acceleration    (gmx_detectcpu_t                data,
                                     FILE *                         log);


#ifdef __cplusplus
}
#endif


#endif /* _GMX_DETECTCPU_H_ */
