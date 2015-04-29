/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
#ifndef GMX_CPUID_H_
#define GMX_CPUID_H_

#include <stdio.h>


#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif


/* Currently identifiable CPU Vendors */
enum gmx_cpuid_vendor
{
    GMX_CPUID_VENDOR_CANNOTDETECT,   /* Should only be used if something fails */
    GMX_CPUID_VENDOR_UNKNOWN,
    GMX_CPUID_VENDOR_INTEL,
    GMX_CPUID_VENDOR_AMD,
    GMX_CPUID_VENDOR_FUJITSU,
    GMX_CPUID_VENDOR_IBM,
    GMX_CPUID_VENDOR_ARM,
    GMX_CPUID_NVENDORS
};


/* CPU feature/property list, to be used as indices into the feature array of the
 * gmxcpuid_t data structure.
 *
 * To facilitate looking things up, we keep this list alphabetical.
 * The list is NOT exhaustive - we have basically added stuff that might be
 * useful in an application like Gromacs.
 *
 * AMD and Intel tend to share most architectural elements, and even if the
 * flags might have to be detected in different ways (different cpuid registers),
 * once the flag is present the functions should be identical. Unfortunately the
 * trend right now (2012) seems to be that they are diverging. This means that
 * we need to use specific flags to the compiler to maximize performance, and
 * then the binaries might not be portable between Intel and AMD as they were
 * before when we only needed to check for SSE and/or SSE2 support in Gromacs.
 */
enum gmx_cpuid_feature
{
    GMX_CPUID_FEATURE_CANNOTDETECT,      /* Flag set if we could not detect on this CPU  */
    GMX_CPUID_FEATURE_X86_AES,           /* x86 advanced encryption standard accel.      */
    GMX_CPUID_FEATURE_X86_APIC,          /* APIC support                                 */
    GMX_CPUID_FEATURE_X86_AVX,           /* Advanced vector extensions                   */
    GMX_CPUID_FEATURE_X86_AVX2,          /* AVX2 including gather support (not used yet) */
    GMX_CPUID_FEATURE_X86_AVX_512F,      /* Foundation AVX-512 instructions              */
    GMX_CPUID_FEATURE_X86_AVX_512PF,     /* Extended gather/scatter for AVX-512          */
    GMX_CPUID_FEATURE_X86_AVX_512ER,     /* Extended-range 1/x and /1sqrt(x) for AVX-512 */
    GMX_CPUID_FEATURE_X86_AVX_512CD,     /* Memory conflict-detection for AVX-512        */
    GMX_CPUID_FEATURE_X86_CLFSH,         /* Supports CLFLUSH instruction                 */
    GMX_CPUID_FEATURE_X86_CMOV,          /* Conditional move insn support                */
    GMX_CPUID_FEATURE_X86_CX8,           /* Supports CMPXCHG8B (8-byte compare-exchange) */
    GMX_CPUID_FEATURE_X86_CX16,          /* Supports CMPXCHG16B (16-byte compare-exchg)  */
    GMX_CPUID_FEATURE_X86_F16C,          /* Supports 16-bit FP conversion instructions   */
    GMX_CPUID_FEATURE_X86_FMA,           /* Fused-multiply add support (mainly for AVX)  */
    GMX_CPUID_FEATURE_X86_FMA4,          /* 4-operand FMA, only on AMD for now           */
    GMX_CPUID_FEATURE_X86_HTT,           /* Hyper-Threading supported                    */
    GMX_CPUID_FEATURE_X86_LAHF_LM,       /* LAHF/SAHF support in 64 bits                 */
    GMX_CPUID_FEATURE_X86_MISALIGNSSE,   /* Support for misaligned SSE data instructions */
    GMX_CPUID_FEATURE_X86_MMX,           /* MMX registers and instructions               */
    GMX_CPUID_FEATURE_X86_MSR,           /* Supports Intel model-specific-registers      */
    GMX_CPUID_FEATURE_X86_NONSTOP_TSC,   /* Invariant TSC (constant rate in ACPI states) */
    GMX_CPUID_FEATURE_X86_PCID,          /* Process context identifier support           */
    GMX_CPUID_FEATURE_X86_PCLMULDQ,      /* Carry-less 64-bit multiplication supported   */
    GMX_CPUID_FEATURE_X86_PDCM,          /* Perfmon and Debug Capability                 */
    GMX_CPUID_FEATURE_X86_PDPE1GB,       /* Support for 1GB pages                        */
    GMX_CPUID_FEATURE_X86_POPCNT,        /* Supports the POPCNT (population count) insn  */
    GMX_CPUID_FEATURE_X86_PSE,           /* Supports 4MB-pages (page size extension)     */
    GMX_CPUID_FEATURE_X86_RDRND,         /* RDRAND high-quality hardware random numbers  */
    GMX_CPUID_FEATURE_X86_RDTSCP,        /* Serializing rdtscp instruction available     */
    GMX_CPUID_FEATURE_X86_SSE2,          /* SSE 2                                        */
    GMX_CPUID_FEATURE_X86_SSE3,          /* SSE 3                                        */
    GMX_CPUID_FEATURE_X86_SSE4A,         /* SSE 4A                                       */
    GMX_CPUID_FEATURE_X86_SSE4_1,        /* SSE 4.1                                      */
    GMX_CPUID_FEATURE_X86_SSE4_2,        /* SSE 4.2                                      */
    GMX_CPUID_FEATURE_X86_SSSE3,         /* Supplemental SSE3                            */
    GMX_CPUID_FEATURE_X86_TDT,           /* TSC deadline timer                           */
    GMX_CPUID_FEATURE_X86_X2APIC,        /* Extended xAPIC Support                       */
    GMX_CPUID_FEATURE_X86_XOP,           /* AMD extended instructions, only AMD for now  */
    GMX_CPUID_FEATURE_ARM_NEON,          /* 32-bit ARM NEON                              */
    GMX_CPUID_FEATURE_ARM_NEON_ASIMD,    /* 64-bit ARM AArch64 Advanced SIMD             */
    GMX_CPUID_FEATURE_IBM_QPX,           /* IBM QPX SIMD (BlueGene/Q and later)          */
    GMX_CPUID_FEATURE_IBM_VMX,           /* IBM VMX SIMD (Altivec on Power6 and later)   */
    GMX_CPUID_FEATURE_IBM_VSX,           /* IBM VSX SIMD (Power7 and later)              */
    GMX_CPUID_NFEATURES
};


/* Currently supported SIMD instruction sets, intrinsics or other similar combinations
 * in Gromacs. There is not always a 1-to-1 correspondence with feature flags; on some AMD
 * hardware we prefer to use 128bit AVX instructions (although 256-bit ones could be executed).
 * These are listed in increasing order for sets supported by one CPU.
 * The order is only used for printing "minimum" and "maximum" suggested
 * SIMD instruction sets for nodes in a cluster, so pairs like
 * GMX_CPUID_SIMD_X86_AVX_128_FMA vs GMX_CPUID_SIMD_X86_AVX_256 which strictly
 * speaking can't be ordered are not really an issue.
 */
enum gmx_cpuid_simd
{
    GMX_CPUID_SIMD_CANNOTDETECT,    /* Should only be used if something fails */
    GMX_CPUID_SIMD_NONE,
    GMX_CPUID_SIMD_REFERENCE,
    GMX_CPUID_SIMD_X86_SSE2,
    GMX_CPUID_SIMD_X86_SSE4_1,
    GMX_CPUID_SIMD_X86_AVX_128_FMA,
    GMX_CPUID_SIMD_X86_AVX_256,
    GMX_CPUID_SIMD_X86_AVX2_256,
    GMX_CPUID_SIMD_X86_AVX_512F,
    GMX_CPUID_SIMD_X86_AVX_512ER,
    GMX_CPUID_SIMD_SPARC64_HPC_ACE,
    GMX_CPUID_SIMD_IBM_QPX,
    GMX_CPUID_SIMD_IBM_VMX,
    GMX_CPUID_SIMD_IBM_VSX,
    GMX_CPUID_SIMD_ARM_NEON,
    GMX_CPUID_SIMD_ARM_NEON_ASIMD,
    GMX_CPUID_NSIMD
};

/* Text strings corresponding to CPU vendors */
extern const char *
gmx_cpuid_vendor_string[GMX_CPUID_NVENDORS];

/* Text strings for CPU feature indices */
extern const char *
gmx_cpuid_feature_string[GMX_CPUID_NFEATURES];

/* Text strings for Gromacs SIMD instruction sets */
extern const char *
gmx_cpuid_simd_string[GMX_CPUID_NSIMD];


/* Abstract data type with CPU detection information. Set by gmx_cpuid_init(). */
typedef struct gmx_cpuid *
    gmx_cpuid_t;


/* Return the SIMD instruction set GROMACS was compiled with. */
enum gmx_cpuid_simd
gmx_compiled_simd           ();


/* Fill the data structure by using CPU detection instructions.
 * Return 0 on success, 1 if something bad happened.
 */
int
gmx_cpuid_init              (gmx_cpuid_t *              cpuid);


/* Return the vendor id as enumerated type. Use gmx_cpuid_vendor_string[]
 * to get the corresponding text string.
 */
enum gmx_cpuid_vendor
gmx_cpuid_vendor            (gmx_cpuid_t                cpuid);


/* Return a constant pointer to the processor brand string. */
const char *
gmx_cpuid_brand             (gmx_cpuid_t                cpuid);


/* Return processor family version. For a chip of version 1.2.3, this is 1 */
int
gmx_cpuid_family            (gmx_cpuid_t                cpuid);

/* Return processor model version, For a chip of version 1.2.3, this is 2. */
int
gmx_cpuid_model             (gmx_cpuid_t                cpuid);

/* Return processor stepping version, For a chip of version 1.2.3, this is 3. */
int
gmx_cpuid_stepping          (gmx_cpuid_t                cpuid);


/* Check whether a particular CPUID feature is set.
 * Returns 0 if flag "feature" is not set, 1 if the flag is set. We cannot use
 * gmx_bool here since this file must be possible to compile without simple.h.
 */
int
gmx_cpuid_feature           (gmx_cpuid_t                cpuid,
                             enum gmx_cpuid_feature     feature);


/* Check whether the CPU is an Intel with Nehalem microarchitecture.
 * Return 0 if not Intel Nehalem, 1 if Intel Nehalem.
 */
int
gmx_cpuid_is_intel_nehalem  (const gmx_cpuid_t          cpuid);


/* Return pointers to cpu topology information.
 *
 * Important: CPU topology requires more OS support than most other
 * functions in this file, including support for thread pinning to hardware.
 * This means it will not work on some platforms, including e.g. Mac OS X.
 * Thus, it is IMPERATIVE that you check the return value from this routine
 * before doing anything with the information. It is only if the return
 * value is zero that the data is valid.
 *
 * For the returned values we have:
 * - nprocessors         Total number of logical processors reported by OS
 * - npackages           Usually number of CPU sockets
 * - ncores_per_package  Number of cores in each package
 * - nhwthreads_per_core Number of hardware threads per core; 2 for hyperthreading.
 * - package_id          Array with the package index for each logical cpu
 * - core_id             Array with local core index for each logical cpu
 * - hwthread_id         Array with local hwthread index for each logical cpu
 * - locality_order      Array with logical cpu numbers, sorted in order
 *                       of physical and logical locality in the system.
 *
 * All arrays are of length nprocessors.
 */
int
gmx_cpuid_topology(gmx_cpuid_t        cpuid,
                   int *              nprocessors,
                   int *              npackages,
                   int *              ncores_per_package,
                   int *              nhwthreads_per_core,
                   const int **       package_id,
                   const int **       core_id,
                   const int **       hwthread_id,
                   const int **       locality_order);

/* Enumerated values for x86 SMT enabled-status. Note that this does not refer
 * to Hyper-Threading support (that is the flag GMX_CPUID_FEATURE_X86_HTT), but
 * whether Hyper-Threading is _enabled_ and _used_ in bios right now.
 */
enum gmx_cpuid_x86_smt
{
    GMX_CPUID_X86_SMT_CANNOTDETECT,
    GMX_CPUID_X86_SMT_DISABLED,
    GMX_CPUID_X86_SMT_ENABLED
};

/* Returns the status of x86 SMT support. IMPORTANT: There are non-zero
 * return values for this routine that still do not indicate supported and
 * enabled smt/Hyper-Threading. You need to carefully check the return value
 * against the enumerated type values to see what you are getting.
 *
 * Long-term, this functionality will move to a new hardware topology detection
 * layer, but that will require a lot of new code and a working interface to the
 * hwloc library. Surprisingly, there is no simple way to find out that
 * Hyper-Threading is actually turned on without fully enumerating and checking
 * all the cores, which we presently can only do on Linux. This means a couple
 * of things:
 *
 * 1) If you want to know whether your CPU _supports_ Hyper-Threading in the
 *    first place, check the GMX_CPUID_FEATURE_X86_HTT flag instead!
 * 2) There are several scenarios where this routine will say that it cannot
 *    detect whether SMT is enabled and used right now.
 * 3) If you need support on non-Linux x86, you have to write it :-)
 * 4) Don't invest too much efforts, since this will be replaced with
 *    full hardware topology detection in the future.
 * 5) Don't worry if the detection does not work. It is not a catastrophe, but
 *    but we get slightly better performance on x86 if we use Hyper-Threading
 *    cores in direct space, but not reciprocal space.
 *
 * Since this routine presently only supports Hyper-Threading we say X86_SMT
 * in order not to give the impression we can detect any SMT. We haven't
 * even tested the performance on other SMT implementations, so it is not
 * obvious we shouldn't use SMT there.
 *
 * Note that you can get more complete topology information from
 * gmx_cpuid_topology(), although that requires slightly more OS support.
 */
enum gmx_cpuid_x86_smt
gmx_cpuid_x86_smt(gmx_cpuid_t cpuid);


/* Formats a text string (up to n characters) from the data structure.
 * The output will have max 80 chars between newline characters.
 */
int
gmx_cpuid_formatstring      (gmx_cpuid_t                cpuid,
                             char *                     s,
                             int                        n);


/* Suggests a suitable gromacs SIMD based on the support in the
 * hardware.
 */
enum gmx_cpuid_simd
gmx_cpuid_simd_suggest  (gmx_cpuid_t                    cpuid);


/* Check if this binary was compiled with the same SIMD instructions as we
 * would suggest for the current hardware. Always print stats to the log file
 * if it is non-NULL, and if we don't have a match, print a warning in log
 * (if non-NULL) and if print_to_stderr!=0 also to stderr.
 * The suggested SIMD instruction set simd_suggest is obtained with
 * gmx_cpuid_simd_suggest(), but with MPI this might be different for
 * different nodes, so it shoul be passed here after parallel reduction.
 */
int
gmx_cpuid_simd_check    (enum gmx_cpuid_simd        simd_suggest,
                         FILE *                     log,
                         int                        print_to_stderr);


/* Release resources used by data structure. Note that the pointer to the
 * CPU brand string will no longer be valid once this routine has been called.
 */
void
gmx_cpuid_done              (gmx_cpuid_t                cpuid);




#ifdef __cplusplus
}
#endif


#endif /* GMX_CPUID_H_ */
