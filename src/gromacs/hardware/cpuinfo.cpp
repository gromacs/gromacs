/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#    include <windows.h>    // sysinfo(), necessary for topology stuff
#endif

#ifdef HAVE_SCHED_H
#    include <sched.h>      // sched_getaffinity(), sched_setaffinity()
#endif
#ifdef HAVE_UNISTD_H
#    include <unistd.h>     // sysconf()
#endif

#include <cctype>
#include <cstdlib>

#include <algorithm>
#include <fstream>
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
void
trimString(std::string * s)
{
    // heading
    s->erase(s->begin(), std::find_if(s->begin(), s->end(), [](char &c) -> bool { return !std::isspace(c); }));
    // trailing
    s->erase(std::find_if(s->rbegin(), s->rend(), [](char &c) -> bool { return !std::isspace(c); }).base(), s->end());
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
int
executeX86CpuID(unsigned int     gmx_unused level,
                unsigned int     gmx_unused ecxval,
                unsigned int *              eax,
                unsigned int *              ebx,
                unsigned int *              ecx,
                unsigned int *              edx)
{
    if (c_architecture == Architecture::X86)
    {
#if defined __GNUC__ || GMX_X86_GCC_INLINE_ASM

        // any compiler that understands gcc inline assembly
        *eax = level;
        *ecx = ecxval;
        *ebx = 0;
        *edx = 0;

#    if (defined __i386__ || defined __i386 || defined _X86_ || defined _M_IX86) && defined(__PIC__)
        // Avoid clobbering the global offset table in 32-bit pic code (ebx register)
        __asm__ __volatile__ ("xchgl %%ebx, %1  \n\t"
                              "cpuid            \n\t"
                              "xchgl %%ebx, %1  \n\t"
                              : "+a" (*eax), "+r" (*ebx), "+c" (*ecx), "+d" (*edx));
#    else
        // i386 without PIC, or x86-64. Things are easy and we can clobber any reg we want
        __asm__ __volatile__ ("cpuid            \n\t"
                              : "+a" (*eax), "+b" (*ebx), "+c" (*ecx), "+d" (*edx));
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

#endif          // check for inline asm on x86
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
 *  If support for the cpuid instruction is present, we check for Intel
 *  or AMD vendors.
 *
 *  \return gmx::CpuInfo::Vendor::Intel, gmx::CpuInfo::Vendor::Amd. If neither
 *          Intel nor Amd can be identified, or if the code fails to execute,
 *          gmx::CpuInfo::Vendor::Unknown is returned.
 */
CpuInfo::Vendor
detectX86Vendor()
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
    }
    return v;
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
void
setFeatureFromBit(std::set<CpuInfo::Feature> *   featureSet,
                  CpuInfo::Feature               feature,
                  unsigned int                   registerValue,
                  unsigned char                  bit)
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
void
detectX86Features(std::string *                  brand,
                  int *                          family,
                  int *                          model,
                  int *                          stepping,
                  std::set<CpuInfo::Feature> *   features)
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

        setFeatureFromBit(features, CpuInfo::Feature::X86_Sse3,     ecx,  0 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Pclmuldq, ecx,  1 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Ssse3,    ecx,  9 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Fma,      ecx, 12 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Cx16,     ecx, 13 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Pdcm,     ecx, 15 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Pcid,     ecx, 17 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Sse4_1,   ecx, 19 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Sse4_2,   ecx, 20 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_X2Apic,   ecx, 21 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Popcnt,   ecx, 23 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Tdt,      ecx, 24 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Aes,      ecx, 25 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx,      ecx, 28 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_F16C,     ecx, 29 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Rdrnd,    ecx, 30 );

        setFeatureFromBit(features, CpuInfo::Feature::X86_Pse,      edx,  3 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Msr,      edx,  5 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Cx8,      edx,  8 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Apic,     edx,  9 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Cmov,     edx, 15 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Clfsh,    edx, 19 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Mmx,      edx, 23 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Sse2,     edx, 26 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Htt,      edx, 28 );
    }

    if (maxStdLevel >= 0x7)
    {
        executeX86CpuID(0x7, 0, &eax, &ebx, &ecx, &edx);

        setFeatureFromBit(features, CpuInfo::Feature::X86_Hle,      ebx,  4 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx2,     ebx,  5 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Rtm,      ebx, 11 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx512F,  ebx, 16 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx512PF, ebx, 26 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx512ER, ebx, 27 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx512CD, ebx, 28 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Sha,      ebx, 29 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx512BW, ebx, 30 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Avx512VL, ebx, 31 );
    }

    // Check whether Hyper-threading is really possible to enable in the hardware,
    // not just technically supported by this generation of processors
    if (features->count(CpuInfo::Feature::X86_Htt) && maxStdLevel >= 0x4)
    {
        executeX86CpuID(0x1, 0, &eax, &ebx, &ecx, &edx);
        unsigned int maxLogicalCores  = (ebx >> 16) & 0x0ff;
        executeX86CpuID(0x4, 0, &eax, &ebx, &ecx, &edx);
        unsigned int maxPhysicalCores = ((eax >> 26) & 0x3f) + 1;
        if (maxLogicalCores/maxPhysicalCores < 2)
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

        setFeatureFromBit(features, CpuInfo::Feature::X86_Lahf,        ecx,  0 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Sse4A,       ecx,  6 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_MisalignSse, ecx,  7 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Xop,         ecx, 11 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Fma4,        ecx, 16 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_PDPE1GB,     edx, 26 );
        setFeatureFromBit(features, CpuInfo::Feature::X86_Rdtscp,      edx, 27 );
    }

    if (maxExtLevel >= 0x80000005)
    {
        // Get the x86 CPU brand string (3 levels, 16 bytes in each)
        brand->clear();
        for (unsigned int level = 0x80000002; level < 0x80000005; level++)
        {
            executeX86CpuID(level, 0, &eax, &ebx, &ecx, &edx);
            // Add eax, ebx, ecx, edx contents as 4 chars each to the brand string
            brand->append(reinterpret_cast<const char *>(&eax), sizeof(eax));
            brand->append(reinterpret_cast<const char *>(&ebx), sizeof(ebx));
            brand->append(reinterpret_cast<const char *>(&ecx), sizeof(ecx));
            brand->append(reinterpret_cast<const char *>(&edx), sizeof(edx));
        }
        trimString(brand);
    }

    if (maxExtLevel >= 0x80000007)
    {
        executeX86CpuID(0x80000007, 0, &eax, &ebx, &ecx, &edx);

        setFeatureFromBit(features, CpuInfo::Feature::X86_NonstopTsc, edx,  8 );
    }
}


/*! \brief Return a vector with x86 APIC IDs for all threads
 *
 *  \param haveX2Apic  True if the processors supports x2APIC, otherwise vanilla APIC.
 *
 *  \returns A new std::vector of unsigned integer APIC IDs, one for each
 *           logical processor in the system.
 */
const std::vector<unsigned int>
detectX86ApicIDs(bool gmx_unused haveX2Apic)
{
    std::vector<unsigned int>  apicID;

    // We cannot just ask for all APIC IDs, but must force execution on each
    // hardware thread and extract the APIC id there.
#if HAVE_SCHED_AFFINITY && defined HAVE_SYSCONF
    unsigned int   eax, ebx, ecx, edx;
    unsigned int   nApic = sysconf(_SC_NPROCESSORS_ONLN);
    cpu_set_t      saveCpuSet;
    cpu_set_t      cpuSet;
    sched_getaffinity(0, sizeof(cpu_set_t), &saveCpuSet);
    CPU_ZERO(&cpuSet);
    for (unsigned int i = 0; i < nApic; i++)
    {
        CPU_SET(i, &cpuSet);
        sched_setaffinity(0, sizeof(cpu_set_t), &cpuSet);
        if (haveX2Apic)
        {
            executeX86CpuID(0xb, 0, &eax, &ebx, &ecx, &edx);
            apicID.push_back(edx);
        }
        else
        {
            executeX86CpuID(0x1, 0, &eax, &ebx, &ecx, &edx);
            apicID.push_back(ebx >> 24);
        }
        CPU_CLR(i, &cpuSet);
    }
    sched_setaffinity(0, sizeof(cpu_set_t), &saveCpuSet);
#elif GMX_NATIVE_WINDOWS
    unsigned int   eax, ebx, ecx, edx;
    SYSTEM_INFO    sysinfo;
    GetSystemInfo( &sysinfo );
    unsigned int   nApic        = sysinfo.dwNumberOfProcessors;
    unsigned int   saveAffinity = SetThreadAffinityMask(GetCurrentThread(), 1);
    for (DWORD_PTR i = 0; i < nApic; i++)
    {
        SetThreadAffinityMask(GetCurrentThread(), (((DWORD_PTR)1)<<i));
        Sleep(0);
        if (haveX2Apic)
        {
            executeX86CpuID(0xb, 0, &eax, &ebx, &ecx, &edx);
            apicID.push_back(edx);
        }
        else
        {
            executeX86CpuID(0x1, 0, &eax, &ebx, &ecx, &edx);
            apicID.push_back(ebx >> 24);
        }
    }
    SetThreadAffinityMask(GetCurrentThread(), saveAffinity);
#endif
    return apicID;
}


/*! \brief Utility to renumber indices extracted from APIC IDs
 *
 * \param v  Vector with unsigned integer indices
 *
 * This routine returns the number of unique different elements found in the vector,
 * and renumbers these starting from 0. For example, the vector {0,1,2,8,9,10,8,9,10,0,1,2}
 * will be rewritten to {0,1,2,3,4,5,3,4,5,0,1,2}, and it returns 6 for the
 * number of unique elements.
 */
void
renumberIndex(std::vector<unsigned int> * v)
{
    std::vector<unsigned int> sortedV (*v);
    std::sort(sortedV.begin(), sortedV.end());

    std::vector<unsigned int> uniqueSortedV (sortedV);
    auto                      it = std::unique(uniqueSortedV.begin(), uniqueSortedV.end());
    uniqueSortedV.resize( std::distance(uniqueSortedV.begin(), it) );

    for (std::size_t i = 0; i < uniqueSortedV.size(); i++)
    {
        unsigned int val = uniqueSortedV[i];
        std::replace_if(v->begin(), v->end(), [val](unsigned int &c) -> bool { return c == val; }, static_cast<unsigned int>(i));
    }
}


/*! \brief Try to detect basic CPU topology information using x86 cpuid
 *
 *  If x2APIC support is present, this is our first choice, otherwise we
 *  attempt to use old vanilla APIC.
 *
 *  \return A new vector of entries with socket, core, hwthread information
 *          for each logical processor.
 */
std::vector<CpuInfo::LogicalProcessor>
detectX86LogicalProcessors()
{
    unsigned int   eax;
    unsigned int   ebx;
    unsigned int   ecx;
    unsigned int   edx;
    unsigned int   maxStdLevel;
    unsigned int   maxExtLevel;
    bool           haveApic;
    bool           haveX2Apic;

    std::vector<CpuInfo::LogicalProcessor> logicalProcessors;

    // Find largest standard & extended level input values allowed
    executeX86CpuID(0x0, 0, &eax, &ebx, &ecx, &edx);
    maxStdLevel = eax;
    executeX86CpuID(0x80000000, 0, &eax, &ebx, &ecx, &edx);
    maxExtLevel = eax;

    if (maxStdLevel >= 0x1)
    {
        executeX86CpuID(0x1, 0, &eax, &ebx, &ecx, &edx);
        haveX2Apic = (ecx & (1 << 21)) && maxStdLevel >= 0xb;
        haveApic   = (edx & (1 <<  9)) && maxExtLevel >= 0x80000008;
    }
    else
    {
        haveX2Apic = false,
        haveApic   = false;
    }

    if (haveX2Apic || haveApic)
    {
        unsigned int   hwThreadBits;
        unsigned int   coreBits;
        // Get bits for cores and hardware threads
        if (haveX2Apic)
        {
            executeX86CpuID(0xb, 0, &eax, &ebx, &ecx, &edx);
            hwThreadBits    = eax & 0x1f;
            executeX86CpuID(0xb, 1, &eax, &ebx, &ecx, &edx);
            coreBits        = (eax & 0x1f) - hwThreadBits;
        }
        else    // haveApic
        {
            // AMD without x2APIC does not support SMT - there are no hwthread bits in apic ID
            hwThreadBits = 0;
            // Get number of core bits in apic ID - try modern extended method first
            executeX86CpuID(0x80000008, 0, &eax, &ebx, &ecx, &edx);
            coreBits = (ecx >> 12) & 0xf;
            if (coreBits == 0)
            {
                // Legacy method for old single/dual core AMD CPUs
                int i = ecx & 0xf;
                while (i >> coreBits)
                {
                    coreBits++;
                }
            }
        }

        std::vector<unsigned int>  apicID = detectX86ApicIDs(haveX2Apic);

        if (!apicID.empty())
        {
            // APIC IDs can be buggy, and it is always a mess. Typically more bits are
            // reserved than needed, and the numbers might not increment by 1 even in
            // a single socket or core. Extract, renumber, and check that things make sense.
            unsigned int               hwThreadMask  = (1 << hwThreadBits) - 1;
            unsigned int               coreMask      = (1 << coreBits) - 1;
            std::vector<unsigned int>  hwThreadRanks;
            std::vector<unsigned int>  coreRanks;
            std::vector<unsigned int>  socketRanks;

            for (auto a : apicID)
            {
                hwThreadRanks.push_back( static_cast<int>( a & hwThreadMask ) );
                coreRanks.push_back( static_cast<int>( ( a >> hwThreadBits ) & coreMask ) );
                socketRanks.push_back( static_cast<int>( a >> ( coreBits + hwThreadBits ) ) );
            }

            renumberIndex(&hwThreadRanks);
            renumberIndex(&coreRanks);
            renumberIndex(&socketRanks);

            unsigned int  hwThreadRankSize = 1 + *std::max_element(hwThreadRanks.begin(), hwThreadRanks.end());
            unsigned int  coreRankSize     = 1 + *std::max_element(coreRanks.begin(), coreRanks.end());
            unsigned int  socketRankSize   = 1 + *std::max_element(socketRanks.begin(), socketRanks.end());

            if (socketRankSize * coreRankSize * hwThreadRankSize == apicID.size() )
            {
                // Alright, everything looks consistent, so put it in the result
                for (std::size_t i = 0; i < apicID.size(); i++)
                {
                    // While the internal APIC IDs are always unsigned integers, we also cast to
                    // plain integers for the externally exposed vectors, since that will make
                    // it possible to use '-1' for invalid entries in the future.
                    logicalProcessors.push_back( { int(socketRanks[i]), int(coreRanks[i]), int(hwThreadRanks[i]) } );
                }
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
const std::map<std::string, std::string>
parseProcCpuInfo()
{
    std::ifstream                       procCpuInfo("/proc/cpuinfo");
    std::string                         line;
    std::map<std::string, std::string>  cpuInfo;

    while (std::getline(procCpuInfo, line))
    {
        if (!line.empty())
        {
            std::stringstream iss(line);
            std::string       key;
            std::string       val;
            std::getline(iss, key, ':');  // part before colon
            std::getline(iss, val);       // part after colon
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
CpuInfo::Vendor
detectProcCpuInfoVendor(const std::map<std::string, std::string> &cpuInfo)
{
    const std::map<std::string, CpuInfo::Vendor> testVendors =
    {
        { "GenuineIntel", CpuInfo::Vendor::Intel   },
        { "Intel",        CpuInfo::Vendor::Intel   },
        { "AuthenticAmd", CpuInfo::Vendor::Amd     },
        { "AMD",          CpuInfo::Vendor::Amd     },
        { "ARM",          CpuInfo::Vendor::Arm     },
        { "AArch64",      CpuInfo::Vendor::Arm     },
        { "Fujitsu",      CpuInfo::Vendor::Fujitsu },
        { "IBM",          CpuInfo::Vendor::Ibm     },
        { "POWER",        CpuInfo::Vendor::Ibm     },
        { "Oracle",       CpuInfo::Vendor::Oracle  },
    };

    // For each label in /proc/cpuinfo, compare the value to the name in the
    // testNames map above, and if it's a match return the vendor.
    for (auto &l : { "vendor_id", "vendor", "manufacture", "model", "processor", "cpu" })
    {
        if (cpuInfo.count(l))
        {
            // there was a line with this left-hand side in /proc/cpuinfo
            const std::string &s1 = cpuInfo.at(l);

            for (auto &t : testVendors)
            {
                const std::string &s2 = t.first;

                // If the entire name we are testing (s2) matches the first part of
                // the string after the colon in /proc/cpuinfo (s1) we found our vendor
                if (std::equal(s2.begin(), s2.end(), s1.begin(),
                               [](const char &x, const char &y) -> bool { return tolower(x) == tolower(y); }))
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
void
detectProcCpuInfoIbm(const std::map<std::string, std::string> &cpuInfo,
                     std::string *                             brand,
                     std::set<CpuInfo::Feature> *              features)
{
    // Get brand string from 'cpu' label if present, otherwise 'Processor'
    if (cpuInfo.count("cpu"))
    {
        *brand = cpuInfo.at("cpu");
    }
    else if (cpuInfo.count("Processor"))
    {
        *brand = cpuInfo.at("Processor");
    }

    if (brand->find("A2") != std::string::npos)
    {
        // If the processor identification contains "A2", this is BlueGene/Q with QPX
        features->insert(CpuInfo::Feature::Ibm_Qpx);
    }

    for (auto &l : { "model name", "model", "Processor", "cpu" })
    {
        if (cpuInfo.count(l))
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
void
detectProcCpuInfoArm(const std::map<std::string, std::string>   &cpuInfo,
                     std::string *                               brand,
                     int *                                       family,
                     int *                                       model,
                     int *                                       stepping,
                     std::set<CpuInfo::Feature> *                features)
{
    if (cpuInfo.count("Processor"))
    {
        *brand = cpuInfo.at("Processor");
    }
    if (cpuInfo.count("CPU architecture"))
    {
        *family = std::strtol(cpuInfo.at("CPU architecture").c_str(), nullptr, 10);
        // For some 64-bit CPUs it appears to say 'AArch64' instead
        if (*family == 0 && cpuInfo.at("CPU architecture").find("AArch64") != std::string::npos)
        {
            *family = 8;  // fragile - no idea how a future ARMv9 will be represented in this case
        }
    }
    if (cpuInfo.count("CPU variant"))
    {
        *model    = std::strtol(cpuInfo.at("CPU variant").c_str(), nullptr, 16);
    }
    if (cpuInfo.count("CPU revision"))
    {
        *stepping = std::strtol(cpuInfo.at("CPU revision").c_str(), nullptr, 10);
    }

    if (cpuInfo.count("Features"))
    {
        const std::string &s = cpuInfo.at("Features");
        if (s.find("neon") != std::string::npos)
        {
            features->insert(CpuInfo::Feature::Arm_Neon);
        }
        if (s.find("asimd") != std::string::npos)
        {
            // At least Jetson TX1 runs a 32-bit environment by default, although
            // the kernel is 64-bits, and reports asimd feature flags. We cannot
            // use Neon-asimd in this case, so make sure we are on a 64-bit platform.
            if (sizeof(void *) == 8)
            {
                features->insert(CpuInfo::Feature::Arm_NeonAsimd);
            }
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
void
detectProcCpuInfo(CpuInfo::Vendor *              vendor,
                  std::string   *                brand,
                  int   *                        family,
                  int   *                        model,
                  int   *                        stepping,
                  std::set<CpuInfo::Feature> *   features)
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
        case CpuInfo::Vendor::Ibm:
            detectProcCpuInfoIbm(cpuInfo, brand, features);
            break;

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
}   // namespace anonymous


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
        detectX86Features(&result.brandString_, &result.family_, &result.model_,
                          &result.stepping_, &result.features_);
        result.logicalProcessors_ = detectX86LogicalProcessors();
    }
    else
    {
        // Not x86
        if (c_architecture == Architecture::Arm)
        {
            result.vendor_  = CpuInfo::Vendor::Arm;
        }
        else if (c_architecture == Architecture::PowerPC)
        {
            result.vendor_  = CpuInfo::Vendor::Ibm;
        }

#if defined __aarch64__ || ( defined _M_ARM && _M_ARM >= 8 )
        result.features_.insert(Feature::Arm_Neon);      // ARMv8 always has Neon
        result.features_.insert(Feature::Arm_NeonAsimd); // ARMv8 always has Neon-asimd
#endif

#if defined sun
        result.vendor_ = CpuInfo::Vendor::Oracle;
#endif

        // On Linux we might be able to find information in /proc/cpuinfo. If vendor or brand
        // is set to a known value this routine will not overwrite it.
        detectProcCpuInfo(&result.vendor_, &result.brandString_, &result.family_,
                          &result.model_, &result.stepping_, &result.features_);
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


CpuInfo::CpuInfo()
    : vendor_(CpuInfo::Vendor::Unknown), brandString_("Unknown CPU brand"),
      family_(0), model_(0), stepping_(0)
{
}


const std::map<CpuInfo::Vendor, std::string>
CpuInfo::s_vendorStrings_ =
{
    { CpuInfo::Vendor::Unknown, "Unknown vendor"                  },
    { CpuInfo::Vendor::Intel, "Intel"                             },
    { CpuInfo::Vendor::Amd, "AMD"                                 },
    { CpuInfo::Vendor::Fujitsu, "Fujitsu"                         },
    { CpuInfo::Vendor::Ibm, "IBM"                                 },
    { CpuInfo::Vendor::Arm, "ARM"                                 },
    { CpuInfo::Vendor::Oracle, "Oracle"                           },
};


const std::map<CpuInfo::Feature, std::string>
CpuInfo::s_featureStrings_ =
{
    { CpuInfo::Feature::X86_Aes, "aes"                            },
    { CpuInfo::Feature::X86_Amd, "amd"                            },
    { CpuInfo::Feature::X86_Apic, "apic"                          },
    { CpuInfo::Feature::X86_Avx, "avx"                            },
    { CpuInfo::Feature::X86_Avx2, "avx2"                          },
    { CpuInfo::Feature::X86_Avx512F, "avx512f"                    },
    { CpuInfo::Feature::X86_Avx512PF, "avx512pf"                  },
    { CpuInfo::Feature::X86_Avx512ER, "avx512er"                  },
    { CpuInfo::Feature::X86_Avx512CD, "avx512cd"                  },
    { CpuInfo::Feature::X86_Avx512BW, "avx512bw"                  },
    { CpuInfo::Feature::X86_Avx512VL, "avx512vl"                  },
    { CpuInfo::Feature::X86_Clfsh, "clfsh"                        },
    { CpuInfo::Feature::X86_Cmov, "cmov"                          },
    { CpuInfo::Feature::X86_Cx8, "cx8"                            },
    { CpuInfo::Feature::X86_Cx16, "cx16"                          },
    { CpuInfo::Feature::X86_F16C, "f16c"                          },
    { CpuInfo::Feature::X86_Fma, "fma"                            },
    { CpuInfo::Feature::X86_Fma4, "fma4"                          },
    { CpuInfo::Feature::X86_Hle, "hle"                            },
    { CpuInfo::Feature::X86_Htt, "htt"                            },
    { CpuInfo::Feature::X86_Intel, "intel"                        },
    { CpuInfo::Feature::X86_Lahf, "lahf"                          },
    { CpuInfo::Feature::X86_MisalignSse, "misalignsse"            },
    { CpuInfo::Feature::X86_Mmx, "mmx"                            },
    { CpuInfo::Feature::X86_Msr, "msr"                            },
    { CpuInfo::Feature::X86_NonstopTsc, "nonstop_tsc"             },
    { CpuInfo::Feature::X86_Pcid, "pcid"                          },
    { CpuInfo::Feature::X86_Pclmuldq, "pclmuldq"                  },
    { CpuInfo::Feature::X86_Pdcm, "pdcm"                          },
    { CpuInfo::Feature::X86_PDPE1GB, "pdpe1gb"                    },
    { CpuInfo::Feature::X86_Popcnt, "popcnt"                      },
    { CpuInfo::Feature::X86_Pse, "pse"                            },
    { CpuInfo::Feature::X86_Rdrnd, "rdrnd"                        },
    { CpuInfo::Feature::X86_Rdtscp, "rdtscp"                      },
    { CpuInfo::Feature::X86_Rtm, "rtm"                            },
    { CpuInfo::Feature::X86_Sha, "sha"                            },
    { CpuInfo::Feature::X86_Sse2, "sse2"                          },
    { CpuInfo::Feature::X86_Sse3, "sse3"                          },
    { CpuInfo::Feature::X86_Sse4A, "sse4a"                        },
    { CpuInfo::Feature::X86_Sse4_1, "sse4.1"                      },
    { CpuInfo::Feature::X86_Sse4_2, "sse4.2"                      },
    { CpuInfo::Feature::X86_Ssse3, "ssse3"                        },
    { CpuInfo::Feature::X86_Tdt, "tdt"                            },
    { CpuInfo::Feature::X86_X2Apic, "x2apic"                      },
    { CpuInfo::Feature::X86_Xop, "xop"                            },
    { CpuInfo::Feature::Arm_Neon, "neon"                          },
    { CpuInfo::Feature::Arm_NeonAsimd, "neon_asimd"               },
    { CpuInfo::Feature::Ibm_Qpx, "qpx"                            },
    { CpuInfo::Feature::Ibm_Vmx, "vmx"                            },
    { CpuInfo::Feature::Ibm_Vsx, "vsx"                            },
    { CpuInfo::Feature::Fujitsu_HpcAce, "hpc-ace"                 }
};


bool
cpuIsX86Nehalem(const CpuInfo &cpuInfo)
{
    return (cpuInfo.vendor() == gmx::CpuInfo::Vendor::Intel &&
            cpuInfo.family() == 6 &&
            (cpuInfo.model() == 0x2E || cpuInfo.model() == 0x1A ||
             cpuInfo.model() == 0x1E || cpuInfo.model() == 0x2F ||
             cpuInfo.model() == 0x2C || cpuInfo.model() == 0x25) );
}

}  // namespace gmx

#ifdef GMX_CPUINFO_STANDALONE
int
main(int argc, char **argv)
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

    std::string   arg(argv[1]);
    gmx::CpuInfo  cpuInfo(gmx::CpuInfo::detect());

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
        for (auto &f : cpuInfo.featureSet() )
        {
            printf("%s ", cpuInfo.featureString(f).c_str());
        }
        printf("\n");
    }
    else if (arg == "-topology")
    {
        // Undocumented debug option, usually not present in standalone version
        for (auto &t : cpuInfo.logicalProcessors() )
        {
            printf("%3u %3u %3u\n", t.socketRankInMachine, t.coreRankInSocket, t.hwThreadRankInCore);
        }
    }
    return 0;
}
#endif
