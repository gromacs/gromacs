/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Declares gmx::CpuInfo
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \inlibraryapi
 * \ingroup module_hardware
 */
#ifndef GMX_HARDWARE_CPUINFO_H
#define GMX_HARDWARE_CPUINFO_H

#include <map>
#include <set>
#include <string>
#include <vector>

namespace gmx
{

/*! \libinternal \brief Detect CPU capabilities and basic logical processor info
 *
 *  This class provides a lot of information about x86 CPUs, and some very
 *  limited information about other hardware. The logical processor information
 *  is only available on x86, and is used as a fallback implementation in
 *  the HardwareTopology class.
 *  If you actually need information about the hardware topology, use the much
 *  more general implementation in the HardwareTopology class instead, since
 *  that will both be more portable and contain more information.
 *
 * \ingroup module_hardware
 */
class CpuInfo
{

public:
    /*! \brief Amount of cpu information present (incremental) */
    enum class SupportLevel
    {
        None,                //!< No cpu information whatsoever. Sorry.
        Name,                //!< Only vendor and/or brand is set
        Features,            //!< Some features are set
        LogicalProcessorInfo //!< Everything includling logical processor information
    };

    /*! \brief Processor/system vendors */
    enum class Vendor
    {
        Unknown,  //!< Unidentified
        Intel,    //!< GenuineIntel
        Amd,      //!< AuthenticAMD
        Fujitsu,  //!< Only works on Linux (parsed from /proc/cpuinfo)
        Ibm,      //!< Only works on Linux (parsed from /proc/cpuinfo)
        Arm,      //!< Only works on Linux (parsed from /proc/cpuinfo)
        Oracle,   //!< Cannot detect anything else yet (no /proc/cpuinfo available)
        Hygon,    //!< HygonGenuine
        RiscV32,  //!< RISC-V 32 bit
        RiscV64,  //!< RISC-V 64 bit
        Loongson, //!< Loongson
    };

    /*! \brief List of CPU features
     *
     *  These values can be used as arguments to the feature() method
     *  to check whether a specific feature was found on the CPU we are
     *  running on.
     */
    enum class Feature
    {
        X86_Aes,             //!< x86 advanced encryption standard accel.
        X86_Amd,             //!< This is an AMD x86 processor
        X86_Apic,            //!< APIC support
        X86_Avx,             //!< Advanced vector extensions
        X86_Avx2,            //!< AVX2 including gather support (not used yet)
        X86_Avx512F,         //!< Foundation AVX-512 instructions
        X86_Avx512PF,        //!< Extended gather/scatter for AVX-512
        X86_Avx512ER,        //!< AVX-512 exponential and reciprocal extensions
        X86_Avx512CD,        //!< Memory conflict-detection for AVX-512
        X86_Avx512BW,        //!< AVX-512 byte and word instructions
        X86_Avx512VL,        //!< AVX-512 vector length extensions
        X86_Avx512BF16,      //!< AVX-512 BFloat16 instructions
        X86_Avx512secondFMA, //!< AVX-512 second FMA unit
        X86_Clfsh,           //!< Supports CLFLUSH instruction
        X86_Cmov,            //!< Conditional move insn support
        X86_Cx8,             //!< Supports CMPXCHG8B (8-byte compare-exchange)
        X86_Cx16,            //!< Supports CMPXCHG16B (16-byte compare-exchg)
        X86_F16C,            //!< Supports 16-bit FP conversion instructions
        X86_Fma,             //!< Fused-multiply add support (mainly for AVX)
        X86_Fma4,            //!< 4-operand FMA, only on AMD for now
        X86_Hle,             //!< Hardware lock elision
        X86_Htt,   //!< Hyper-Threading enabled (NOTE: might not match the CPUID HTT support flag)
        X86_Intel, //!< This is an Intel x86 processor
        X86_Lahf,  //!< LAHF/SAHF support in 64 bits
        X86_MisalignSse, //!< Support for misaligned SSE data instructions
        X86_Mmx,         //!< MMX registers and instructions
        X86_Msr,         //!< Supports Intel model-specific-registers
        X86_NonstopTsc,  //!< Invariant TSC (constant rate in ACPI states)
        X86_Pcid,        //!< Process context identifier support
        X86_Pclmuldq,    //!< Carry-less 64-bit multiplication supported
        X86_Pdcm,        //!< Perfmon and Debug Capability
        X86_PDPE1GB,     //!< Support for 1GB pages
        X86_Popcnt,      //!< Supports the POPCNT (population count) insn
        X86_Pse,         //!< Supports 4MB-pages (page size extension)
        X86_Rdrnd,       //!< RDRAND high-quality hardware random numbers
        X86_Rdtscp,      //!< Serializing rdtscp instruction available
        X86_Rtm,         //!< Restricted transactional memory
        X86_Sha,         //!< Intel SHA extensions
        X86_Sse2,        //!< SSE 2
        X86_Sse3,        //!< SSE 3
        X86_Sse4A,       //!< SSE 4A
        X86_Sse4_1,      //!< SSE 4.1
        X86_Sse4_2,      //!< SSE 4.2
        X86_Ssse3,       //!< Supplemental SSE3
        X86_Tdt,         //!< TSC deadline timer
        X86_X2Apic,      //!< Extended xAPIC Support
        X86_Xop,         //!< AMD extended instructions, only AMD for now
        Arm_Neon,        //!< 32-bit ARM NEON
        Arm_NeonAsimd,   //!< 64-bit ARM AArch64 Advanced SIMD
        Arm_Sve,         //!< ARM Scalable Vector Extensions
        Ibm_Qpx,         //!< IBM QPX SIMD (BlueGene/Q)
        Ibm_Vmx,         //!< IBM VMX SIMD (Altivec on Power6 and later)
        Ibm_Vsx,         //!< IBM VSX SIMD (Power7 and later)
        Fujitsu_HpcAce,  //!< Fujitsu Sparc64 HPC-ACE
        X86_Hygon        //!< This is a Hygon x86 processor
    };

    /*! \libinternal \brief Entry with basic information for a single processing unit.
     *
     * To avoid duplicating functionality with the hardware topology, be aware
     * that the ids in this structure correspond to the raw APIC (or similar) ids
     * reported by the hardware. These are not necessarily ranks, so some
     * ids might not be present, and others might be larger than the total
     * number of such elements present - be careful when using them.
     */
    struct LogicalProcessor
    {
        int packageIdInMachine; //!< Id of the current package in the system
        int coreIdInPackage;    //!< Id of the current core in its package
        int puIdInCore;         //!< Id of logical cpu inside its core
        int osId;               //!< Id corresponding to OS logical cpu enumeration
    };

    /*! \brief Perform detection and construct a CpuInfo class from the results.
     *
     *  \note The detection should generally be performed again in different
     *        contexts.  This might seem like overkill, but there
     *        are systems (e.g. Arm) where processors can go completely offline
     *        during deep sleep, so at least in theory it is good to have a
     *        possibility of forcing re-detection if necessary.
     */
    static CpuInfo detect();

    /*! \brief Check what cpu information is available
     *
     *  The amount of cpu information that can be detected depends on the
     *  OS, compiler, and CPU, and on non-x86 platforms it can be fragile.
     *  Before basing decisions on the output or warning the user about
     *  optimizations, you want to check whether it was possible to detect
     *  the information you need.
     */
    SupportLevel supportLevel() const { return supportLevel_; }

    /*! \brief Enumerated value for vendor */
    Vendor vendor() const { return vendor_; }

    /*! \brief String description of vendor:
     *
     *  \throws std::out_of_range if the vendor is not present in the internal
     *          map of vendor names. This can only happen if we extend the enum
     *          type but forget to add the string with the vendor name.
     */
    const std::string& vendorString() const;

    /*! \brief String description of processor */
    const std::string& brandString() const { return brandString_; }

    /*! \brief Major version/generation of the processor */
    int family() const { return family_; }

    /*! \brief Middle version of the processor */
    int model() const { return model_; }

    /*! \brief Minor version of the processor */
    int stepping() const { return stepping_; }

    /*! \brief Check for availability of specific feature
     *
     *  \param f  feature to query support for
     *
     *  \return True if the feature is available, otherwise false.
     */
    bool feature(Feature f) const
    {
        // If the entry is present in the set it is supported
        return (features_.count(f) != 0);
    }

    /*! \brief String description of a specific feature
     *
     *  \throws std::out_of_range if the feature is not present in the internal
     *          map of feature names. This can only happen if we extend the enum
     *          type but forget to add the string with the feature name.
     */
    static const std::string& featureString(Feature f);

    /*! \brief Set of all supported features on this processor
     *
     *  This is only intended for logfiles, debugging or similar output when we
     *  need a full list of all the features available on the CPU.
     */
    const std::set<Feature>& featureSet() const { return features_; }

    /*! \brief Reference to processing unit topology
     *
     *  Only a few systems (e.g. x86) provide logical processor information in cpuinfo.
     *  This method returns a reference to a vector, whose length will either be
     *  zero (if topology information is not available) or the number of
     *  processing units on which we were allowed to execute, as defined
     *  e.g. by cpu sets.
     *
     *  This is only meant to be use as a fallback implementation for our
     *  HardwareTopology class; any user code that needs access to hardware
     *  topology information should use that class instead.
     *
     *  \note For clarity, it is likely better to use the supportLevel()
     *        method to check if this information is available rather than
     *        relying on the length of the vector.
     */
    const std::vector<LogicalProcessor>& logicalProcessors() const { return logicalProcessors_; }

private:
    CpuInfo();

    SupportLevel                  supportLevel_;      //!< Available cpuinfo information
    Vendor                        vendor_;            //!<  Value of vendor for current cpu
    std::string                   brandString_;       //!<  Text description of cpu
    int                           family_;            //!<  Major version of current cpu
    int                           model_;             //!<  Middle version of current cpu
    int                           stepping_;          //!<  Minor version of current cpu
    std::set<Feature>             features_;          //!< Set of features supported on this cpu
    std::vector<LogicalProcessor> logicalProcessors_; //!< Simple logical processor topology
};                                                    // class CpuInfo

/*! \brief Return true if the CPU is an Intel x86 Nehalem
 *
 * \param cpuInfo  Object with cpu information
 *
 * \returns  True if running on Nehalem CPU
 */
bool cpuIsX86Nehalem(const CpuInfo& cpuInfo);

/*! \brief Return true if the CPU is a first generation AMD Zen (produced by AMD or Hygon)
 *
 * \param cpuInfo  Object with cpu information
 *
 * \returns  True if running on a first generation AMD Zen
 */
bool cpuIsAmdZen1(const CpuInfo& cpuInfo);

} // namespace gmx

#endif // GMX_HARDWARE_CPUINFO_H
