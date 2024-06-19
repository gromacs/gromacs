/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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

/*! \internal \file
 * \brief
 * Implements gmx::CpuInfo.
 *
 * We need to be able to compile this file in stand-alone mode to use basic
 * CPU feature detection to set the SIMD acceleration and similar things in
 * CMake, while we still want to use more features that enable topology
 * detection when config.h is present.
 *
 * We solve this by skipping the advanced stuff when the preprocessor
 * macro GMX_CPUINFO_STANDALONE is defined. In this case you likely also need to
 * define GMX_X86_GCC_INLINE_ASM if you are on x86; without inline assembly
 * support it is not possible to perform the actual detection on Linux/Mac.
 * Since these macros are specific to this file, they do not use the GMX prefix.
 *
 * The remaining defines (GMX_NATIVE_WINDOWS,HAVE_UNISTD_H,HAVE_SCHED_H,
 * HAVE_SYSCONF, HAVE_SCHED_AFFINITY) are only used to determine the topology on
 * 86, and for this we rely on including config.h.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_hardware
 */

#ifndef GMX_CPUINFO_STANDALONE
#    include "gmxpre.h"
#endif

#include "cpuinfo.h"

#ifndef GMX_CPUINFO_STANDALONE
#    include "config.h"
#else
#    define GMX_NATIVE_WINDOWS 0
#endif

#if defined _MSC_VER
#    include <intrin.h> // __cpuid()
#endif

#if GMX_NATIVE_WINDOWS
#    include <windows.h> // sysinfo(), necessary for topology stuff
#endif

#ifdef HAVE_SCHED_H
#    include <sched.h> // sched_getaffinity(), sched_setaffinity()
#endif
#ifdef HAVE_UNISTD_H
#    include <unistd.h> // sysconf()
#endif

#include <cctype>
#include <cstdint> // uint32_t in X86 processor name code
#include <cstdlib>

#include <algorithm>
#include <fstream>
#include <initializer_list>
#include <map>
#include <set>
#include <sstream>
#include <string>

#ifdef GMX_CPUINFO_STANDALONE
#    define gmx_unused
#else
#    include "gromacs/utility/basedefinitions.h"
#endif

#include "architecture.h"

//! Convenience macro to help us avoid ifdefs each time we use sysconf
#if !defined(_SC_NPROCESSORS_ONLN) && defined(_SC_NPROC_ONLN)
#    define _SC_NPROCESSORS_ONLN _SC_NPROC_ONLN
#endif

namespace gmx
{

namespace
{

/*! \cond internal */

/******************************************************************************
 *                                                                            *
 *   Utility functions to make this file independent of the GROMACS library   *
 *                                                                            *
 ******************************************************************************/

/*! \brief Remove initial and trailing whitespace from string
 *
 *  \param s  Pointer to string where whitespace will be removed
 */
void trimString(std::string* s)
{
    // heading
    s->erase(s->begin(),
             std::find_if(s->begin(), s->end(), [](char& c) -> bool { return std::isspace(c) == 0; }));
    // trailing
    s->erase(
            std::find_if(s->rbegin(), s->rend(), [](char& c) -> bool { return std::isspace(c) == 0; })
                    .base(),
            s->end());
}


/******************************************************************************
 *                                                                            *
 *                         x86 detection functions                            *
 *                                                                            *
 ******************************************************************************/

/*! \brief execute x86 cpuid instructions with custom level and extended level
 *
 *  \param level   The main cpuid level (input argument for eax register)
 *  \param ecxval  Extended level (input argument for ecx register)
 *  \param eax     Output in eax register
 *  \param ebx     Output in ebx register
 *  \param ecx     Output in ecx register
 *  \param edx     Output in edx register
 *
 *  \return 0 on success, or non-zero if the instruction could not execute.
 */
int executeX86CpuID(unsigned int gmx_unused level,
                    unsigned int gmx_unused ecxval,
                    unsigned int*           eax,
                    unsigned int*           ebx,
                    unsigned int*           ecx,
                    unsigned int*           edx)
{
    if (c_architecture == Architecture::X86)
    {
#if defined __GNUC__ || GMX_X86_GCC_INLINE_ASM

        // any compiler that understands gcc inline assembly
        *eax = level;
        *ecx = ecxval;
        *ebx = 0;
        *edx = 0;

#    if GMX_IS_X86_32 && defined(__PIC__)
        // Avoid clobbering the global offset table in 32-bit pic code (ebx register)
        __asm__ __volatile__(
                "xchgl %%ebx, %1  \n\t"
                "cpuid            \n\t"
                "xchgl %%ebx, %1  \n\t"
                : "+a"(*eax), "+r"(*ebx), "+c"(*ecx), "+d"(*edx));
#    elif GMX_IS_X86_64
        // i386 without PIC, or x86-64. Things are easy and we can clobber any reg we want
        __asm__ __volatile__("cpuid            \n\t"
                             : "+a"(*eax), "+b"(*ebx), "+c"(*ecx), "+d"(*edx));
#    else
        // Not a normal x86, which could happen when a compiler
        // targetting non-x86 pretends to be GCC.
#    endif
        return 0;

#elif defined _MSC_VER

        // MSVC (and icc on windows) on ia32 or x86-64
        int cpuInfo[4];
        __cpuidex(cpuInfo, level, ecxval);
        *eax = static_cast<unsigned int>(cpuInfo[0]);
        *ebx = static_cast<unsigned int>(cpuInfo[1]);
        *ecx = static_cast<unsigned int>(cpuInfo[2]);
        *edx = static_cast<unsigned int>(cpuInfo[3]);
        return 0;

#else

        // We are on x86, but without compiler support for cpuid if we get here
        *eax = 0;
        *ebx = 0;
        *ecx = 0;
        *edx = 0;
        return 1;

#endif // check for inline asm on x86
    }
    else
    {
        // We are not on x86
        *eax = 0;
        *ebx = 0;
        *ecx = 0;
        *edx = 0;
        return 1;
    }
}


/*! \brief Detect x86 vendors by using the cpuid assembly instructions
 *
 *  If support for the cpuid instruction is present, we check for Intel,
 *  AMD or Hygon vendors
 *
 *  \return gmx::CpuInfo::Vendor::Intel, gmx::CpuInfo::Vendor::Amd,
 *          gmx::CpuInfl::Vendor::Hygon, . If neither Intel, Amd  nor
 *          Hygon can be identified, or if the code fails to execute,
 *          gmx::CpuInfo::Vendor::Unknown is returned.
 */
CpuInfo::Vendor detectX86Vendor()
{
    unsigned int    eax, ebx, ecx, edx;
    CpuInfo::Vendor v = CpuInfo::Vendor::Unknown;

    if (executeX86CpuID(0x0, 0, &eax, &ebx, &ecx, &edx) == 0)
    {
        if (ebx == 0x756e6547 && ecx == 0x6c65746e && edx == 0x49656e69)
        {
            v = CpuInfo::Vendor::Intel; // ebx=='uneG', ecx=='letn', edx=='Ieni'
        }
        else if (ebx == 0x68747541 && ecx == 0x444d4163 && edx == 0x69746e65)
        {
            v = CpuInfo::Vendor::Amd; // ebx=='htuA', ecx=='DMAc', edx=='itne'
        }
        else if (ebx == 0x6f677948 && ecx == 0x656e6975 && edx == 0x6e65476e)
        {
            v = CpuInfo::Vendor::Hygon; // ebx=='ogyH', ecx=='eniu', edx=='neGn'
        }
    }
    return v;
}

/*! \brief Detect second AVX-512 FMA from the processor name
 *
 * Should only be called for processors already determined to support AVX-512.
 *
 *  \param [in] brand     x86 processor name
 *  \param [in] model     x86 model
 *  \return               True if second FMA present
 */
bool detectProcCpuInfoSecondAvx512FMA(const std::string& brand, int model)
{
    // Skylake server
    if (model == 0x55)
    {
        // detect Xeon
        if (brand.find("Xeon") == 9)
        {
            // detect Silver or Bronze or specific models
            if (brand.find("Silver") == 17 || brand.find("Bronze") == 17
                || (brand.find('W') == 17 && brand.find('0') == 21)   // detect Xeon W 210x
                || (brand.find('D') == 17 && brand.find("21") == 19)) // detect Xeon D 2xxx
            {
                return false;
            }
            // detect Gold 5xxx - can be corrected once Cooper Lake is added
            else if (brand.find("Gold") == 17 && brand.find('5') == 22)
            {
                return (brand.find("53") == 22 || // detect Cooper Lake
                        brand.find("22") == 24);  // detect 5[12]22
            }
        }
        return true;
    }
    // Cannon Lake client
    if (model == 0x66)
    {
        return false;
    }
    // Ice Lake client
    if (model == 0x7d || model == 0x7e)
    {
        return false;
    }
    // This is the right default...
    return true;
}

/*! \brief Simple utility function to set/clear feature in a set
 *
 *  \param featureSet    Pointer to the feature set to update
 *  \param feature       The specific feature to set/clear
 *  \param registerValue Register value (returned from cpuid)
 *  \param bit           Bit to check in registerValue. The feature will be
 *                       added to the featureSet if this bit is set.
 *
 *  \note Nothing is done if the bit is not set. In particular, this will not
 *        erase anything if the feature already exists in the set.
 */
void setFeatureFromBit(std::set<CpuInfo::Feature>* featureSet,
                       CpuInfo::Feature            feature,
                       unsigned int                registerValue,
                       unsigned char               bit)
{
    if (registerValue & (1 << bit))
    {
        featureSet->insert(feature);
    }
}

/*! \brief Process x86 cpuinfo features that are common to Intel and AMD CPUs
 *
 *  \param[out] brand      String where to write the x86 brand string
 *  \param[out] family     Major version of processor
 *  \param[out] model      Middle version of processor
 *  \param[out] stepping   Minor version of processor
 *  \param[out] features   Feature set where supported features are inserted
 */
void detectX86Features(std::string* brand, int* family, int* model, int* stepping, std::set<CpuInfo::Feature>* features)
{
    unsigned int eax, ebx, ecx, edx;

    // Return if we cannot execute any levels
    if (executeX86CpuID(0x0, 0, &eax, &ebx, &ecx, &edx) != 0)
    {
        return;
    }
    unsigned int maxStdLevel = eax;

    if (maxStdLevel >= 0x1)
    {
        executeX86CpuID(0x1, 0, &eax, &ebx, &ecx, &edx);

        *family   = ((eax & 0x0ff00000) >> 20) + ((eax & 0x00000f00) >> 8);
        *model    = ((eax & 0x000f0000) >> 12) + ((eax & 0x000000f0) >> 4);
        *stepping = (eax & 0x0000000f);

        setFeatureFromBit(features, CpuInfo::Feature::X86_Sse3, ecx, 0);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Pclmuldq, ecx, 1);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Ssse3, ecx, 9);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Fma, ecx, 12);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Cx16, ecx, 13);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Pdcm, ecx, 15);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Pcid, ecx, 17);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Sse4_1, ecx, 19);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Sse4_2, ecx, 20);
        setFeatureFromBit(features, CpuInfo::Feature::X86_X2Apic, ecx, 21);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Popcnt, ecx, 23);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Tdt, ecx, 24);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Aes, ecx, 25);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx, ecx, 28);
        setFeatureFromBit(features, CpuInfo::Feature::X86_F16C, ecx, 29);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Rdrnd, ecx, 30);

        setFeatureFromBit(features, CpuInfo::Feature::X86_Pse, edx, 3);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Msr, edx, 5);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Cx8, edx, 8);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Apic, edx, 9);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Cmov, edx, 15);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Clfsh, edx, 19);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Mmx, edx, 23);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Sse2, edx, 26);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Htt, edx, 28);
    }

    // Check whether Hyper-threading is really possible to enable in the hardware,
    // not just technically supported by this generation of processors
    if ((features->count(CpuInfo::Feature::X86_Htt) != 0U) && maxStdLevel >= 0x4)
    {
        executeX86CpuID(0x1, 0, &eax, &ebx, &ecx, &edx);
        unsigned int maxLogicalCores = (ebx >> 16) & 0x0ff;
        executeX86CpuID(0x4, 0, &eax, &ebx, &ecx, &edx);
        unsigned int maxPhysicalCores = ((eax >> 26) & 0x3f) + 1;
        if (maxLogicalCores / maxPhysicalCores < 2)
        {
            features->erase(CpuInfo::Feature::X86_Htt);
        }
    }

    if (executeX86CpuID(0x80000000, 0, &eax, &ebx, &ecx, &edx) != 0)
    {
        // No point in continuing if we don't support any extended levels
        return;
    }
    unsigned int maxExtLevel = eax;

    if (maxExtLevel >= 0x80000001)
    {
        executeX86CpuID(0x80000001, 0, &eax, &ebx, &ecx, &edx);

        setFeatureFromBit(features, CpuInfo::Feature::X86_Lahf, ecx, 0);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Sse4A, ecx, 6);
        setFeatureFromBit(features, CpuInfo::Feature::X86_MisalignSse, ecx, 7);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Xop, ecx, 11);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Fma4, ecx, 16);
        setFeatureFromBit(features, CpuInfo::Feature::X86_PDPE1GB, edx, 26);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Rdtscp, edx, 27);
    }

    if (maxExtLevel >= 0x80000005)
    {
        // Get the x86 CPU brand string (3 levels, 16 bytes in each)
        brand->clear();
        for (unsigned int level = 0x80000002; level < 0x80000005; level++)
        {
            executeX86CpuID(level, 0, &eax, &ebx, &ecx, &edx);
            // Add eax, ebx, ecx, edx contents as 4 chars each to the brand string
            brand->append(reinterpret_cast<const char*>(&eax), sizeof(eax));
            brand->append(reinterpret_cast<const char*>(&ebx), sizeof(ebx));
            brand->append(reinterpret_cast<const char*>(&ecx), sizeof(ecx));
            brand->append(reinterpret_cast<const char*>(&edx), sizeof(edx));
        }
        trimString(brand);
    }

    if (maxStdLevel >= 0x7)
    {
        executeX86CpuID(0x7, 0, &eax, &ebx, &ecx, &edx);

        setFeatureFromBit(features, CpuInfo::Feature::X86_Hle, ebx, 4);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx2, ebx, 5);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Rtm, ebx, 11);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx512F, ebx, 16);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx512PF, ebx, 26);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx512ER, ebx, 27);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx512CD, ebx, 28);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Sha, ebx, 29);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx512BW, ebx, 30);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx512VL, ebx, 31);

        executeX86CpuID(0x7, 0x1, &eax, &ebx, &ecx, &edx);
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx512BF16, eax, 5);

        if (features->count(CpuInfo::Feature::X86_Avx512F) != 0)
        {
            // Only checking if the CPU supports AVX-512. There is no CPUID bit for this.
            if (detectProcCpuInfoSecondAvx512FMA(*brand, *model))
            {
                features->insert(CpuInfo::Feature::X86_Avx512secondFMA);
            }
        }
    }


    if (maxExtLevel >= 0x80000007)
    {
        executeX86CpuID(0x80000007, 0, &eax, &ebx, &ecx, &edx);

        setFeatureFromBit(features, CpuInfo::Feature::X86_NonstopTsc, edx, 8);
    }
}

/*! \brief internal structure to return os logical cpu id together with APIC info for it */
struct ApicInfo
{
    int          osId;   //!< OS index of logical cpu/processor
    unsigned int apicId; //!< APIC id obtained from cpuid when running on the osId processor
};

/*! \brief Return a vector with x86 APIC info for all processing units
 *
 *  \param haveX2Apic  True if the processors supports x2APIC, otherwise vanilla APIC.
 *
 *  \returns A new vector with os-provided logical cpu id and corresponding hardware APIC id.
 */
std::vector<ApicInfo> detectX86ApicInfo(bool gmx_unused haveX2Apic)
{
    std::vector<ApicInfo> apicInfo;

    // We cannot just ask for all APIC IDs, but must force execution on each
    // hardware thread and extract the APIC id there.
    // Since the thread might have been pinned, we cannot use the present
    // CPU set, but must try all possible CPUs and check the return code
    // if we were allowed to run on that CPU, and eventually restore
    // the original cpu set.
#if HAVE_SCHED_AFFINITY && defined HAVE_SYSCONF
    unsigned int eax, ebx, ecx, edx;
    unsigned int nApic = sysconf(_SC_NPROCESSORS_ONLN);
    cpu_set_t    saveCpuSet;
    sched_getaffinity(0, sizeof(cpu_set_t), &saveCpuSet);
    // We only test threads up to the number of CPUs in the system, not the
    // (much larger) max value of the data structures.
    for (unsigned int i = 0; i < nApic; i++)
    {
        cpu_set_t cpuSet;
        CPU_ZERO(&cpuSet);
        CPU_SET(i, &cpuSet);
        if (sched_setaffinity(0, sizeof(cpu_set_t), &cpuSet) == 0)
        {
            // 0 means the affinity could be set
            if (haveX2Apic)
            {
                executeX86CpuID(0xb, 0, &eax, &ebx, &ecx, &edx);
                apicInfo.push_back({ static_cast<int>(i), edx });
            }
            else
            {
                executeX86CpuID(0x1, 0, &eax, &ebx, &ecx, &edx);
                apicInfo.push_back({ static_cast<int>(i), ebx >> 24 });
            }
        }
    }
    sched_setaffinity(0, sizeof(cpu_set_t), &saveCpuSet);
#elif GMX_NATIVE_WINDOWS
    unsigned int eax, ebx, ecx, edx;
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);

    // calling SetThreadAffinityMask returns the current affinity mask if it succeeds,
    // otherwise zero - so try all processors until we are successful with one so
    // we can save the original affinity mask.
    unsigned int saveThreadAffinity = 0;
    for (DWORD_PTR i = 0; i < sysinfo.dwNumberOfProcessors && saveThreadAffinity == 0; i++)
    {
        saveThreadAffinity = SetThreadAffinityMask(GetCurrentThread(), (((DWORD_PTR)1) << i));
    }

    for (DWORD_PTR i = 0; i < sysinfo.dwNumberOfProcessors; i++)
    {
        if (SetThreadAffinityMask(GetCurrentThread(), (((DWORD_PTR)1) << i)) != 0)
        {
            // On windows, the call will have returned the old mask (non-zero) if it succeeded
            if (haveX2Apic)
            {
                executeX86CpuID(0xb, 0, &eax, &ebx, &ecx, &edx);
                apicInfo.push_back({ static_cast<int>(i), edx });
            }
            else
            {
                executeX86CpuID(0x1, 0, &eax, &ebx, &ecx, &edx);
                apicInfo.push_back({ static_cast<int>(i), ebx >> 24 });
            }
        }
    }
    SetThreadAffinityMask(GetCurrentThread(), saveThreadAffinity);
#endif
    return apicInfo;
}

/*! \brief The layout of the bits in the APIC ID */
struct ApicIdLayout
{
    unsigned int hwThreadBits; //!< The number of least significant bits for hw-threads
    unsigned int coreBits;     //!< The number of core bits following the  hw-thread bits
};

/*! \brief Detect the APIC ID layout for x2APIC
 */
ApicIdLayout detectX2ApicIdLayout()
{
    ApicIdLayout layout;

    unsigned int eax;
    unsigned int ebx;
    unsigned int ecx;
    unsigned int edx;
    executeX86CpuID(0xb, 0, &eax, &ebx, &ecx, &edx);
    layout.hwThreadBits = eax & 0x1f;
    executeX86CpuID(0xb, 1, &eax, &ebx, &ecx, &edx);
    layout.coreBits = (eax & 0x1f) - layout.hwThreadBits;

    return layout;
}

/*! \brief Detect the APIC ID layout for standard APIC or xAPIC on AMD
 *
 * \param[in] maxExtLevel  The largest CPUID extended function input value supported by the processor implementation
 */
ApicIdLayout detectAmdApicIdLayout(unsigned int maxExtLevel)
{
    ApicIdLayout layout;

    unsigned int eax;
    unsigned int ebx;
    unsigned int ecx;
    unsigned int edx;
    executeX86CpuID(0x1, 0, &eax, &ebx, &ecx, &edx);
    int family = ((eax & 0x0ff00000) >> 20) + ((eax & 0x00000f00) >> 8);
    executeX86CpuID(0x80000001, 0, &eax, &ebx, &ecx, &edx);
    bool haveExtendedTopology = (ecx & (1 << 22)) != 0U;

    // NOTE: Here we assume 1 thread per core, unless we have family >= 17h
    layout.hwThreadBits = 0;
    if (family >= 0x17 && haveExtendedTopology && maxExtLevel >= 0x8000001e)
    {
        executeX86CpuID(0x8000001e, 1, &eax, &ebx, &ecx, &edx);
        int numThreadsPerCore = ((ebx >> 8) & 0xff) + 1;
        // NOTE: The AMD documentation only specifies the layout of apicid
        //       when we have 1 or 2 threads per core.
        while (numThreadsPerCore > (1 << layout.hwThreadBits))
        {
            layout.hwThreadBits++;
        }
    }

    // Get number of core bits in apic ID - try modern extended method first
    executeX86CpuID(0x80000008, 0, &eax, &ebx, &ecx, &edx);
    layout.coreBits = (ecx >> 12) & 0xf;
    if (layout.coreBits == 0)
    {
        // Legacy method for old single/dual core AMD CPUs
        int i = ecx & 0xf;
        while (i >> layout.coreBits)
        {
            layout.coreBits++;
        }
    }

    return layout;
}

/*! \brief Try to detect basic CPU topology information using x86 cpuid
 *
 *  If x2APIC support is present, this is our first choice, otherwise we
 *  attempt to use old vanilla APIC.
 *
 *  \return A new vector of entries with package, core, processing unit information
 *          for each logical processor.
 */
std::vector<CpuInfo::LogicalProcessor> detectX86LogicalProcessors()
{
    unsigned int eax;
    unsigned int ebx;
    unsigned int ecx;
    unsigned int edx;
    unsigned int maxStdLevel;
    unsigned int maxExtLevel;
    bool         haveApic;
    bool         haveX2Apic;

    std::vector<CpuInfo::LogicalProcessor> logicalProcessors;

    // Find largest standard & extended level input values allowed
    executeX86CpuID(0x0, 0, &eax, &ebx, &ecx, &edx);
    maxStdLevel = eax;
    executeX86CpuID(0x80000000, 0, &eax, &ebx, &ecx, &edx);
    maxExtLevel = eax;

    if (maxStdLevel >= 0x1)
    {
        executeX86CpuID(0x1, 0, &eax, &ebx, &ecx, &edx);
        haveX2Apic = ((ecx & (1 << 21)) != 0U) && maxStdLevel >= 0xb;
        haveApic   = ((edx & (1 << 9)) != 0U) && maxExtLevel >= 0x80000008;
    }
    else
    {
        haveX2Apic = false;
        haveApic   = false;
    }

    if (haveX2Apic || haveApic)
    {
        ApicIdLayout layout;
        // Get bits for cores and hardware threads
        if (haveX2Apic)
        {
            layout = detectX2ApicIdLayout();
        }
        else // haveApic
        {
            if (detectX86Vendor() == CpuInfo::Vendor::Amd || detectX86Vendor() == CpuInfo::Vendor::Hygon)
            {
                layout = detectAmdApicIdLayout(maxExtLevel);

                if (layout.hwThreadBits > 1)
                {
                    // At the time of writing this code we do not know what
                    // to do with more than 2 threads, so return empty.
                    return logicalProcessors;
                }
            }
            else
            {
                // We do not know the APIC ID layout, return empty.
                return logicalProcessors;
            }
        }

        std::vector<ApicInfo> apicInfo = detectX86ApicInfo(haveX2Apic);
        if (!apicInfo.empty())
        {
            // APIC IDs can be buggy, and it is always a mess. Typically more bits are
            // reserved than needed, and the numbers might not increment by 1 even in
            // a single package or core.
            // We sheepishly avoid dealing with that here, but just extract the information
            // and leave the renumbering to the hardware topology.
            unsigned int hwThreadMask = (1 << layout.hwThreadBits) - 1;
            unsigned int coreMask     = (1 << layout.coreBits) - 1;

            for (const auto& a : apicInfo)
            {
                // Note that these labels are arbitrary and local to the package/core;
                // they will in general NOT correspond to what we later have in the hardware topology
                int packageId = a.apicId >> (layout.coreBits + layout.hwThreadBits);
                int coreId    = (a.apicId >> layout.hwThreadBits) & coreMask;
                int puId      = a.apicId & hwThreadMask;
                logicalProcessors.push_back({ packageId, coreId, puId, a.osId });
            }
        }
    }
    return logicalProcessors; // Will only have contents if everything worked
}


/******************************************************************************
 *                                                                            *
 *              Generic Linux detection by parsing /proc/cpuinfo              *
 *                                                                            *
 ******************************************************************************/

/*! \brief Parse /proc/cpuinfo into a simple string map
 *
 * This routine will read the contents of /proc/cpuinfo, and for each
 * line that is not empty we will assign the (trimmed) string to the right of
 * the colon as a key, and the left-hand side as the value in the map.
 * For multi-processor systems where lines are repeated the latter lines will
 * overwrite the first occurrence.
 *
 * \return New map with the contents. If the file is not available, the returned
 *         map will be empty.
 */
std::map<std::string, std::string> parseProcCpuInfo()
{
    std::ifstream                      procCpuInfo("/proc/cpuinfo");
    std::string                        line;
    std::map<std::string, std::string> cpuInfo;

    while (std::getline(procCpuInfo, line))
    {
        if (!line.empty())
        {
            std::stringstream iss(line);
            std::string       key;
            std::string       val;
            std::getline(iss, key, ':'); // part before colon
            std::getline(iss, val);      // part after colon
            trimString(&key);
            trimString(&val);
            // put it in the map. This will overwrite previous processors, but we don't care.
            cpuInfo[key] = val;
        }
    }
    return cpuInfo;
}


/*! \brief Try to detect vendor from /proc/cpuinfo
 *
 *  \param cpuInfo  Map returned from parseProcCpuinfo()
 *
 *  This routine tries to match a few common labels in /proc/cpuinfo to see if
 *  they begin with the name of a standard vendor. If the file cannot be read
 *  or if no match is found, we return gmx::CpuInfo::Vendor::Unknown.
 */
CpuInfo::Vendor detectProcCpuInfoVendor(const std::map<std::string, std::string>& cpuInfo)
{
    const std::map<std::string, CpuInfo::Vendor> testVendors = {
        { "GenuineIntel", CpuInfo::Vendor::Intel },
        { "Intel", CpuInfo::Vendor::Intel },
        { "AuthenticAmd", CpuInfo::Vendor::Amd },
        { "AMD", CpuInfo::Vendor::Amd },
        { "ARM", CpuInfo::Vendor::Arm },
        { "AArch64", CpuInfo::Vendor::Arm },
        { "Fujitsu", CpuInfo::Vendor::Fujitsu },
        { "IBM", CpuInfo::Vendor::Ibm },
        { "POWER", CpuInfo::Vendor::Ibm },
        { "Oracle", CpuInfo::Vendor::Oracle },
        { "HygonGenuine", CpuInfo::Vendor::Hygon },
        { "Hygon", CpuInfo::Vendor::Hygon },
        { "riscv64", CpuInfo::Vendor::RiscV64 },
        { "riscv32", CpuInfo::Vendor::RiscV32 },
        { "riscv", CpuInfo::Vendor::RiscV32 }, // Must come after riscv64 to avoid misidentification
        { "Loongson", CpuInfo::Vendor::Loongson },
        { "loongarch64", CpuInfo::Vendor::Loongson },
        { "loong64", CpuInfo::Vendor::Loongson },
    };

    // For each label in /proc/cpuinfo, compare the value to the name in the
    // testNames map above, and if it's a match return the vendor.
    for (const auto& l :
         { "vendor_id", "vendor", "manufacture", "model", "processor", "cpu", "Architecture" })
    {
        if (cpuInfo.count(l) != 0U)
        {
            // there was a line with this left-hand side in /proc/cpuinfo
            const std::string& s1 = cpuInfo.at(l);

            for (const auto& t : testVendors)
            {
                const std::string& s2 = t.first;

                // If the entire name we are testing (s2) matches the first part of
                // the string after the colon in /proc/cpuinfo (s1) we found our vendor
                if (std::equal(s2.begin(), s2.end(), s1.begin(), [](const char& x, const char& y) -> bool {
                        return tolower(x) == tolower(y);
                    }))
                {
                    return t.second;
                }
            }
        }
    }
    return CpuInfo::Vendor::Unknown;
}


/*! \brief Detect IBM processor name and features from /proc/cpuinfo
 *
 *  \param      cpuInfo    Map returned from parseProcCpuinfo()
 *  \param[out] brand      String where to write the brand string
 *  \param[out] features   Feature set where supported features are inserted
 *
 *  This routine tries to match a few common labels in /proc/cpuinfo to see if
 *  we can find the processor name and features. It is likely fragile.
 */
void detectProcCpuInfoIbm(const std::map<std::string, std::string>& cpuInfo,
                          std::string*                              brand,
                          std::set<CpuInfo::Feature>*               features)
{
    // Get brand string from 'cpu' label if present, otherwise 'Processor'
    if (cpuInfo.count("cpu") != 0U)
    {
        *brand = cpuInfo.at("cpu");
    }
    else if (cpuInfo.count("Processor") != 0U)
    {
        *brand = cpuInfo.at("Processor");
    }

    if (brand->find("A2") != std::string::npos)
    {
        // If the processor identification contains "A2", this is BlueGene/Q with QPX
        features->insert(CpuInfo::Feature::Ibm_Qpx);
    }

    for (const auto& l : { "model name", "model", "Processor", "cpu" })
    {
        if (cpuInfo.count(l) != 0U)
        {
            std::string s1 = cpuInfo.at(l);
            std::transform(s1.begin(), s1.end(), s1.begin(), ::tolower);

            if (s1.find("altivec") != std::string::npos)
            {
                features->insert(CpuInfo::Feature::Ibm_Vmx);
                // If this is a power6, we only have VMX. All later processors have VSX.
                if (s1.find("power6") == std::string::npos)
                {
                    features->insert(CpuInfo::Feature::Ibm_Vsx);
                }
            }
        }
    }
}


/*! \brief Detect ARM processor name and features from /proc/cpuinfo
 *
 *  \param      cpuInfo    Map returned from parseProcCpuinfo()
 *  \param[out] brand      String where to write the brand string
 *  \param[out] family     Major version of processor
 *  \param[out] model      Middle version of processor
 *  \param[out] stepping   Minor version of processor
 *  \param[out] features   Feature set where supported features are inserted
 *
 *  This routine tries to match a few common labels in /proc/cpuinfo to see if
 *  we can find the processor name and features. It is likely fragile.
 */
void detectProcCpuInfoArm(const std::map<std::string, std::string>& cpuInfo,
                          std::string*                              brand,
                          int*                                      family,
                          int*                                      model,
                          int*                                      stepping,
                          std::set<CpuInfo::Feature>*               features)
{
    if (cpuInfo.count("Processor") != 0U)
    {
        *brand = cpuInfo.at("Processor");
    }
    else if (cpuInfo.count("model name") != 0U)
    {
        *brand = cpuInfo.at("model name");
    }
    else if (cpuInfo.count("CPU part") != 0U)
    {
        // If all else fails, at least print the hexadecimal CPU part code.
        // Mapping can be found at https://github.com/torvalds/linux/blob/master/arch/arm64/include/asm/cputype.h
        *brand = cpuInfo.at("CPU part");
    }

    if (cpuInfo.count("CPU architecture") != 0U)
    {
        *family = std::strtol(cpuInfo.at("CPU architecture").c_str(), nullptr, 10);
        // For some 64-bit CPUs it appears to say 'AArch64' instead
        if (*family == 0 && cpuInfo.at("CPU architecture").find("AArch64") != std::string::npos)
        {
            *family = 8; // fragile - no idea how a future ARMv9 will be represented in this case
        }
    }
    if (cpuInfo.count("CPU variant") != 0U)
    {
        *model = std::strtol(cpuInfo.at("CPU variant").c_str(), nullptr, 16);
    }
    if (cpuInfo.count("CPU revision") != 0U)
    {
        *stepping = std::strtol(cpuInfo.at("CPU revision").c_str(), nullptr, 10);
    }

    if (cpuInfo.count("Features") != 0U)
    {
        const std::string& s = cpuInfo.at("Features");
        if (s.find("neon") != std::string::npos)
        {
            features->insert(CpuInfo::Feature::Arm_Neon);
        }
        if (s.find("asimd") != std::string::npos)
        {
            // At least Jetson TX1 runs a 32-bit environment by default, although
            // the kernel is 64-bits, and reports asimd feature flags. We cannot
            // use Neon-asimd in this case, so make sure we are on a 64-bit platform.
            if (sizeof(void*) == 8)
            {
                features->insert(CpuInfo::Feature::Arm_NeonAsimd);
            }
        }
        if (s.find("sve") != std::string::npos)
        {
            features->insert(CpuInfo::Feature::Arm_Sve);
        }
    }
}


/*! \brief Try to detect vendor, cpu and features from /proc/cpuinfo
 *
 *  \param[out] vendor     Detected hardware vendor
 *  \param[out] brand      String where to write the brand string
 *  \param[out] family     Major version of processor
 *  \param[out] model      Middle version of processor
 *  \param[out] stepping   Minor version of processor
 *  \param[out] features   Feature set where supported features are inserted
 *
 *  This routine reads the /proc/cpuinfo file into a map and calls subroutines
 *  that attempt to parse by matching keys and values to known strings. It is
 *  much more fragile than our x86 detection, but it does not depend on
 *  specific system calls, intrinsics or assembly instructions.
 */
void detectProcCpuInfo(CpuInfo::Vendor*            vendor,
                       std::string*                brand,
                       int*                        family,
                       int*                        model,
                       int*                        stepping,
                       std::set<CpuInfo::Feature>* features)
{
    std::map<std::string, std::string> cpuInfo = parseProcCpuInfo();

    if (*vendor == CpuInfo::Vendor::Unknown)
    {
        *vendor = detectProcCpuInfoVendor(cpuInfo);
    }

    // Unfortunately there is no standard for contents in /proc/cpuinfo. We cannot
    // indiscriminately look for e.g. 'cpu' since it could be either name or an index.
    // To handle this slightly better we use one subroutine per vendor.
    switch (*vendor)
    {
        case CpuInfo::Vendor::Ibm: detectProcCpuInfoIbm(cpuInfo, brand, features); break;

        case CpuInfo::Vendor::Arm:
            detectProcCpuInfoArm(cpuInfo, brand, family, model, stepping, features);
            break;

        default:
            // We only have a single check for fujitsu for now
#ifdef __HPC_ACE__
            features->insert(CpuInfo::Feature::Fujitsu_HpcAce);
#endif
            break;
    }
}
/*! \endcond */
} // namespace


// static
CpuInfo CpuInfo::detect()
{
    CpuInfo result;

    if (c_architecture == Architecture::X86)
    {
        result.vendor_ = detectX86Vendor();

        if (result.vendor_ == CpuInfo::Vendor::Intel)
        {
            result.features_.insert(CpuInfo::Feature::X86_Intel);
        }
        else if (result.vendor_ == CpuInfo::Vendor::Amd)
        {
            result.features_.insert(CpuInfo::Feature::X86_Amd);
        }
        else if (result.vendor_ == CpuInfo::Vendor::Hygon)
        {
            result.features_.insert(CpuInfo::Feature::X86_Hygon);
        }
        detectX86Features(
                &result.brandString_, &result.family_, &result.model_, &result.stepping_, &result.features_);
        result.logicalProcessors_ = detectX86LogicalProcessors();
    }
    else
    {
        // Not x86
        if (c_architecture == Architecture::Arm)
        {
            result.vendor_ = CpuInfo::Vendor::Arm;
        }
        else if (c_architecture == Architecture::PowerPC)
        {
            result.vendor_ = CpuInfo::Vendor::Ibm;
        }
        else if (c_architecture == Architecture::RiscV32)
        {
            result.vendor_ = CpuInfo::Vendor::RiscV32;
        }
        else if (c_architecture == Architecture::RiscV64)
        {
            result.vendor_ = CpuInfo::Vendor::RiscV64;
        }
        else if (c_architecture == Architecture::Loongarch64)
        {
            result.vendor_ = CpuInfo::Vendor::Loongson;
        }

#if defined __aarch64__ || (defined _M_ARM && _M_ARM >= 8)
        result.features_.insert(Feature::Arm_Neon);      // ARMv8 always has Neon
        result.features_.insert(Feature::Arm_NeonAsimd); // ARMv8 always has Neon-asimd
#endif
#if defined __arch64__ && defined __ARM_FEATURE_SVE
        result.features_.insert(Feature::Arm_Sve);
#endif

#if defined sun
        result.vendor_ = CpuInfo::Vendor::Oracle;
#endif

        // On Linux we might be able to find information in /proc/cpuinfo. If vendor or brand
        // is set to a known value this routine will not overwrite it.
        detectProcCpuInfo(&result.vendor_,
                          &result.brandString_,
                          &result.family_,
                          &result.model_,
                          &result.stepping_,
                          &result.features_);
    }

    if (!result.logicalProcessors_.empty())
    {
        result.supportLevel_ = CpuInfo::SupportLevel::LogicalProcessorInfo;
    }
    else if (!result.features_.empty())
    {
        result.supportLevel_ = CpuInfo::SupportLevel::Features;
    }
    else if (result.vendor_ != CpuInfo::Vendor::Unknown
             || result.brandString_ != "Unknown CPU brand")
    {
        result.supportLevel_ = CpuInfo::SupportLevel::Name;
    }
    else
    {
        result.supportLevel_ = CpuInfo::SupportLevel::None;
    }

    return result;
}

CpuInfo::CpuInfo() :
    vendor_(CpuInfo::Vendor::Unknown), brandString_("Unknown CPU brand"), family_(0), model_(0), stepping_(0)
{
}

const std::string& CpuInfo::vendorString() const
{
    static const std::map<Vendor, std::string> vendorStrings = {
        { Vendor::Unknown, "Unknown vendor" },
        { Vendor::Intel, "Intel" },
        { Vendor::Amd, "AMD" },
        { Vendor::Fujitsu, "Fujitsu" },
        { Vendor::Ibm, "IBM" },
        { Vendor::Arm, "ARM" },
        { Vendor::Oracle, "Oracle" },
        { Vendor::Hygon, "Hygon" },
        { Vendor::RiscV32, "RISC-V 32" },
        { Vendor::RiscV64, "RISC-V 64" },
        { Vendor::Loongson, "Loongson" },
    };

    return vendorStrings.at(vendor_);
}


const std::string& CpuInfo::featureString(Feature f)
{
    static const std::map<Feature, std::string> featureStrings = {
        { Feature::X86_Aes, "aes" },
        { Feature::X86_Amd, "amd" },
        { Feature::X86_Apic, "apic" },
        { Feature::X86_Avx, "avx" },
        { Feature::X86_Avx2, "avx2" },
        { Feature::X86_Avx512F, "avx512f" },
        { Feature::X86_Avx512PF, "avx512pf" },
        { Feature::X86_Avx512ER, "avx512er" },
        { Feature::X86_Avx512CD, "avx512cd" },
        { Feature::X86_Avx512BW, "avx512bw" },
        { Feature::X86_Avx512VL, "avx512vl" },
        { Feature::X86_Avx512BF16, "avx512bf16" },
        { Feature::X86_Avx512secondFMA, "avx512secondFMA" },
        { Feature::X86_Clfsh, "clfsh" },
        { Feature::X86_Cmov, "cmov" },
        { Feature::X86_Cx8, "cx8" },
        { Feature::X86_Cx16, "cx16" },
        { Feature::X86_F16C, "f16c" },
        { Feature::X86_Fma, "fma" },
        { Feature::X86_Fma4, "fma4" },
        { Feature::X86_Hle, "hle" },
        { Feature::X86_Htt, "htt" },
        { Feature::X86_Intel, "intel" },
        { Feature::X86_Lahf, "lahf" },
        { Feature::X86_MisalignSse, "misalignsse" },
        { Feature::X86_Mmx, "mmx" },
        { Feature::X86_Msr, "msr" },
        { Feature::X86_NonstopTsc, "nonstop_tsc" },
        { Feature::X86_Pcid, "pcid" },
        { Feature::X86_Pclmuldq, "pclmuldq" },
        { Feature::X86_Pdcm, "pdcm" },
        { Feature::X86_PDPE1GB, "pdpe1gb" },
        { Feature::X86_Popcnt, "popcnt" },
        { Feature::X86_Pse, "pse" },
        { Feature::X86_Rdrnd, "rdrnd" },
        { Feature::X86_Rdtscp, "rdtscp" },
        { Feature::X86_Rtm, "rtm" },
        { Feature::X86_Sha, "sha" },
        { Feature::X86_Sse2, "sse2" },
        { Feature::X86_Sse3, "sse3" },
        { Feature::X86_Sse4A, "sse4a" },
        { Feature::X86_Sse4_1, "sse4.1" },
        { Feature::X86_Sse4_2, "sse4.2" },
        { Feature::X86_Ssse3, "ssse3" },
        { Feature::X86_Tdt, "tdt" },
        { Feature::X86_X2Apic, "x2apic" },
        { Feature::X86_Xop, "xop" },
        { Feature::Arm_Neon, "neon" },
        { Feature::Arm_NeonAsimd, "neon_asimd" },
        { Feature::Arm_Sve, "sve" },
        { Feature::Ibm_Qpx, "qpx" },
        { Feature::Ibm_Vmx, "vmx" },
        { Feature::Ibm_Vsx, "vsx" },
        { Feature::Fujitsu_HpcAce, "hpc-ace" },
        { Feature::X86_Hygon, "hygon" }
    };
    return featureStrings.at(f);
}


bool cpuIsX86Nehalem(const CpuInfo& cpuInfo)
{
    return (cpuInfo.vendor() == CpuInfo::Vendor::Intel && cpuInfo.family() == 6
            && (cpuInfo.model() == 0x2E || cpuInfo.model() == 0x1A || cpuInfo.model() == 0x1E
                || cpuInfo.model() == 0x2F || cpuInfo.model() == 0x2C || cpuInfo.model() == 0x25));
}

bool cpuIsAmdZen1(const CpuInfo& cpuInfo)
{
    /* Both Zen/Zen+/Zen2 have family==23
     * Model numbers for Zen:
     * 1)  Naples, Whitehaven, Summit Ridge, and Snowy Owl;
     * 17) Raven Ridge.
     * Model numbers for Zen+:
     * 8)  Pinnacle Ridge;
     * 24) Picasso.
     * Hygon got license for Zen1, but not Zen2 (https://www.tomshardware.com/news/amd-zen-china-x86-ip-license,39573.html)
     */
    return (cpuInfo.vendor() == CpuInfo::Vendor::Amd && cpuInfo.family() == 23
            && (cpuInfo.model() == 1 || cpuInfo.model() == 17 || cpuInfo.model() == 8
                || cpuInfo.model() == 24))
           || (cpuInfo.vendor() == CpuInfo::Vendor::Hygon);
}

} // namespace gmx

#ifdef GMX_CPUINFO_STANDALONE
int main(int argc, char** argv)
{
    if (argc < 2)
    {
        fprintf(stdout,
                "Usage:\n\n%s [flags]\n\n"
                "Available flags:\n"
                "-vendor        Print CPU vendor.\n"
                "-brand         Print CPU brand string.\n"
                "-family        Print CPU family version.\n"
                "-model         Print CPU model version.\n"
                "-stepping      Print CPU stepping version.\n"
                "-features      Print CPU feature flags.\n",
                argv[0]);
        exit(1);
    }

    std::string  arg(argv[1]);
    gmx::CpuInfo cpuInfo(gmx::CpuInfo::detect());

    if (arg == "-vendor")
    {
        printf("%s\n", cpuInfo.vendorString().c_str());
    }
    else if (arg == "-brand")
    {
        printf("%s\n", cpuInfo.brandString().c_str());
    }
    else if (arg == "-family")
    {
        printf("%d\n", cpuInfo.family());
    }
    else if (arg == "-model")
    {
        printf("%d\n", cpuInfo.model());
    }
    else if (arg == "-stepping")
    {
        printf("%d\n", cpuInfo.stepping());
    }
    else if (arg == "-features")
    {
        // Separate the feature strings with spaces. Note that in the
        // GROMACS cmake code, surrounding whitespace is first
        // stripped by the CPU detection routine, and then added back
        // in the code for making the SIMD suggestion.
        for (const auto& f : cpuInfo.featureSet())
        {
            printf("%s ", cpuInfo.featureString(f).c_str());
        }
        printf("\n");
    }
    else if (arg == "-topology")
    {
        // Undocumented debug option, usually not present in standalone version
        printf("// logical processors we were allowed to run on, mapped to\n"
               "// packageIdInMachine, coreIdInPackage, and puIdInCore\n"
               "{\n");
        for (const auto& t : cpuInfo.logicalProcessors())
        {
            printf("    { %3d , { %3d, %3d, %3d } },\n", t.osId, t.packageIdInMachine, t.coreIdInPackage, t.puIdInCore);
        }
        printf("};\n");
    }
    return 0;
}
#endif
