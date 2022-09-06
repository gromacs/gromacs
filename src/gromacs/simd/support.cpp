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

/*! \internal \file
 *
 * \brief Implements SIMD architecture support query routines
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \ingroup module_simd
 */

#include "gmxpre.h"

#include "gromacs/simd/support.h"

#include "config.h"

#if GMX_SIMD_ARM_SVE
#    include <arm_sve.h>
#endif

#include <cstdio>
#include <cstdlib>

#include <map>
#include <string>

#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/hardware/identifyavx512fmaunits.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/*! \cond libapi */

const std::string& simdString(SimdType s)
{
    static const std::map<SimdType, std::string> name = { { SimdType::None, "None" },
                                                          { SimdType::Reference, "Reference" },
                                                          { SimdType::Generic, "Generic" },
                                                          { SimdType::X86_Sse2, "SSE2" },
                                                          { SimdType::X86_Sse4_1, "SSE4.1" },
                                                          { SimdType::X86_Avx128Fma, "AVX_128_FMA" },
                                                          { SimdType::X86_Avx, "AVX_256" },
                                                          { SimdType::X86_Avx2, "AVX2_256" },
                                                          { SimdType::X86_Avx2_128, "AVX2_128" },
                                                          { SimdType::X86_Avx512, "AVX_512" },
                                                          { SimdType::X86_Avx512Knl, "AVX_512_KNL" },
                                                          { SimdType::Arm_NeonAsimd,
                                                            "ARM_NEON_ASIMD" },
                                                          { SimdType::Arm_Sve, "ARM_SVE" },
                                                          { SimdType::Ibm_Vsx, "IBM_VSX" } };

    return name.at(s);
}

SimdType simdSuggested(const CpuInfo& c)
{
    SimdType suggested = SimdType::None;

    if (c.supportLevel() >= CpuInfo::SupportLevel::Features)
    {
        switch (c.vendor())
        {
            case CpuInfo::Vendor::Intel:
                if (c.feature(CpuInfo::Feature::X86_Avx512ER))
                {
                    suggested = SimdType::X86_Avx512Knl;
                }
                else if (c.feature(CpuInfo::Feature::X86_Avx512F))
                {
                    // If we could not identify the number of AVX512 FMA units we assume 2
                    suggested = (identifyAvx512FmaUnits() == 1) ? SimdType::X86_Avx2 : SimdType::X86_Avx512;
                }
                else if (c.feature(CpuInfo::Feature::X86_Avx2))
                {
                    suggested = SimdType::X86_Avx2;
                }
                else if (c.feature(CpuInfo::Feature::X86_Avx))
                {
                    suggested = SimdType::X86_Avx;
                }
                else if (c.feature(CpuInfo::Feature::X86_Sse4_1))
                {
                    suggested = SimdType::X86_Sse4_1;
                }
                else if (c.feature(CpuInfo::Feature::X86_Sse2))
                {
                    suggested = SimdType::X86_Sse2;
                }
                break;
            case CpuInfo::Vendor::Amd:
            case CpuInfo::Vendor::Hygon:
                if (c.feature(CpuInfo::Feature::X86_Avx2))
                {
                    // AMD Zen supports 256-bit AVX2, but Zen1 performs better with 128-bit
                    // since it can execute two independent such instructions per cycle,
                    // and wider SIMD has slightly lower efficiency in GROMACS.
                    // However... Zen2 supports full-width execution of 256-bit AVX2,
                    // so we only want to apply this hack to Zen/Zen+.
                    suggested = cpuIsAmdZen1(c) ? SimdType::X86_Avx2_128 : SimdType::X86_Avx2;
                }
                else if (c.feature(CpuInfo::Feature::X86_Avx))
                {
                    // Use 128-bit FMA4 SIMD if Fma4 flag is set, otherwise plain 256-bit AVX
                    if (c.feature(CpuInfo::Feature::X86_Fma4))
                    {
                        suggested = SimdType::X86_Avx128Fma;
                    }
                    else
                    {
                        suggested = SimdType::X86_Avx;
                    }
                }
                else if (c.feature(CpuInfo::Feature::X86_Sse4_1))
                {
                    suggested = SimdType::X86_Sse4_1;
                }
                else if (c.feature(CpuInfo::Feature::X86_Sse2))
                {
                    suggested = SimdType::X86_Sse2;
                }

                break;
            case CpuInfo::Vendor::Arm:
                if (c.feature(CpuInfo::Feature::Arm_Sve))
                {
                    suggested = SimdType::Arm_Sve;
                }
                else if (c.feature(CpuInfo::Feature::Arm_NeonAsimd))
                {
                    suggested = SimdType::Arm_NeonAsimd;
                }
                else if (c.feature(CpuInfo::Feature::Arm_Neon))
                {
                    suggested = SimdType::None;
                }
                break;
            case CpuInfo::Vendor::Ibm:
                if (c.feature(CpuInfo::Feature::Ibm_Vsx))
                {
                    suggested = SimdType::Ibm_Vsx;
                }
                else if (c.feature(CpuInfo::Feature::Ibm_Vmx))
                {
                    suggested = SimdType::None;
                }
                break;
            case CpuInfo::Vendor::Fujitsu:
                if (c.feature(CpuInfo::Feature::Fujitsu_HpcAce))
                {
                    suggested = SimdType::None;
                }
                break;
            default: break;
        }
    }
    return suggested;
}

SimdType simdCompiled()
{
#if GMX_SIMD_X86_AVX_512_KNL
    return SimdType::X86_Avx512Knl;
#elif GMX_SIMD_X86_AVX_512
    return SimdType::X86_Avx512;
#elif GMX_SIMD_X86_AVX2_256
    return SimdType::X86_Avx2;
#elif GMX_SIMD_X86_AVX2_128
    return SimdType::X86_Avx2_128;
#elif GMX_SIMD_X86_AVX_256
    return SimdType::X86_Avx;
#elif GMX_SIMD_X86_AVX_128_FMA
    return SimdType::X86_Avx128Fma;
#elif GMX_SIMD_X86_SSE4_1
    return SimdType::X86_Sse4_1;
#elif GMX_SIMD_X86_SSE2
    return SimdType::X86_Sse2;
#elif GMX_SIMD_ARM_NEON_ASIMD
    return SimdType::Arm_NeonAsimd;
#elif GMX_SIMD_ARM_SVE
    return SimdType::Arm_Sve;
#elif GMX_SIMD_IBM_VSX
    return SimdType::Ibm_Vsx;
#elif GMX_SIMD_REFERENCE
    return SimdType::Reference;
#else
    return SimdType::None;
#endif
}

bool simdCheck(gmx::SimdType wanted, FILE* log, bool warnToStdErr)
{
    SimdType compiled = simdCompiled();

    gmx::TextLineWrapper wrapper;
    std::string          logMsg;
    std::string          warnMsg;

    wrapper.settings().setLineLength(78);

    if (compiled == SimdType::X86_Avx2 && wanted == SimdType::X86_Avx512)
    {
        logMsg  = wrapper.wrapToString(formatString(
                "Highest SIMD level supported by all nodes in run: %s\n"
                "SIMD instructions selected at compile time:       %s\n"
                "This program was compiled for different hardware than you are running on, "
                "which could influence performance. This build might have been configured on "
                "a login node with only a single AVX-512 FMA unit (in which case AVX2 is faster), "
                "while the node you are running on has dual AVX-512 FMA units.",
                simdString(wanted).c_str(),
                simdString(compiled).c_str()));
        warnMsg = wrapper.wrapToString(formatString(
                "Compiled SIMD: %s, but for this host/run %s might be better (see log).",
                simdString(compiled).c_str(),
                simdString(wanted).c_str()));
    }
    else if (compiled == SimdType::X86_Avx512 && wanted == SimdType::X86_Avx2
             && identifyAvx512FmaUnits() == 1)
    {
        // The reason for explicitly checking the number of FMA units above is to avoid triggering
        // this conditional if the AVX2 SIMD was supported by some other node in a heterogeneous MPI run.
        logMsg  = wrapper.wrapToString(formatString(
                "Highest SIMD level supported by all nodes in run: %s\n"
                "SIMD instructions selected at compile time:       %s\n"
                "This program was compiled for different hardware than you are running on, "
                "which could influence performance."
                "This host supports AVX-512, but since it only has 1 AVX-512"
                "FMA unit, it would be faster to use AVX2 instead.",
                simdString(wanted).c_str(),
                simdString(compiled).c_str()));
        warnMsg = wrapper.wrapToString(formatString(
                "Compiled SIMD: %s, but for this host/run %s might be better (see log).",
                simdString(compiled).c_str(),
                simdString(wanted).c_str()));
    }
    else if (compiled == SimdType::X86_Avx2 && wanted == SimdType::X86_Avx2_128)
    {
        // Wanted SimdType::X86_Avx2_128 can only be the AMD Zen architecture.
        // AVX2_256 is only up to a few percent slower than AVX2_128
        // in both single and double precision. AVX2_256 is slightly
        // faster with nonbondeds and PME on a GPU. Don't warn the user.
    }
    else if (compiled > wanted && !(compiled == SimdType::X86_Avx && wanted == SimdType::X86_Avx128Fma)
             && !(compiled == SimdType::X86_Avx2_128 && wanted == SimdType::X86_Avx2))
    {
        // Normally it is close to catastrophic if the compiled SIMD type is larger than
        // the supported one, but AVX128Fma is an exception: AMD CPUs will (strongly) prefer
        // AVX128Fma, but they will work fine with AVX too. Thus, make an exception for this.
        // Similarly, AVX2_256 and AVX2_128 are the same instruction set.
        logMsg = wrapper.wrapToString(
                formatString("Highest SIMD level supported by all nodes in run: %s\n"
                             "SIMD instructions selected at compile time:       %s\n"
                             "Compiled SIMD likely not supported by hardware; program might crash.",
                             simdString(wanted).c_str(),
                             simdString(compiled).c_str()));
        warnMsg = logMsg;
    }
    else if (compiled == SimdType::X86_Avx128Fma && wanted != compiled)
    {
        logMsg  = wrapper.wrapToString(formatString(
                "Highest SIMD level supported by all nodes in run: %s\n"
                "SIMD instructions selected at compile time:       %s\n"
                "AMD's early FMA4 AVX extension do not work on modern CPUs; program might crash.",
                simdString(wanted).c_str(),
                simdString(compiled).c_str()));
        warnMsg = logMsg;
    }
    else if (wanted != compiled)
    {
        // This warning will also occur if compiled is X86_Avx and wanted is X86_Avx128Fma
        logMsg  = wrapper.wrapToString(formatString(
                "Highest SIMD level supported by all nodes in run: %s\n"
                "SIMD instructions selected at compile time:       %s\n"
                "This program was compiled for different hardware than you are running on, "
                "which could influence performance.",
                simdString(wanted).c_str(),
                simdString(compiled).c_str()));
        warnMsg = wrapper.wrapToString(formatString(
                "Compiled SIMD: %s, but for this host/run %s might be better (see log).",
                simdString(compiled).c_str(),
                simdString(wanted).c_str()));
#if GMX_SIMD_ARM_SVE
    }
    else if ((compiled == SimdType::Arm_Sve) && (svcntb() != GMX_SIMD_ARM_SVE_LENGTH_VALUE / 8))
    {
        logMsg  = wrapper.wrapToString(formatString(
                "Longest SVE length supported by all nodes in run: %d\n"
                "SVE length selected at compile time:               %ld\n"
                "This program was compiled for different hardware than you are running on, "
                "which will lead to incorrect behavior.\n"
                "Aborting",
                GMX_SIMD_ARM_SVE_LENGTH_VALUE,
                svcntb() * 8));
        warnMsg = wrapper.wrapToString(formatString(
                "Compiled SVE Length: %d, but for this process requires %ld (see log).",
                GMX_SIMD_ARM_SVE_LENGTH_VALUE,
                svcntb() * 8));
#endif
    }

    if (!logMsg.empty() && log != nullptr)
    {
        fprintf(log, "%s\n", logMsg.c_str());
    }
    if (!warnMsg.empty() && warnToStdErr)
    {
        fprintf(stderr, "%s\n", warnMsg.c_str());
    }
#if GMX_SIMD_ARM_SVE
    if ((compiled == SimdType::Arm_Sve) && (svcntb() != GMX_SIMD_ARM_SVE_LENGTH_VALUE / 8))
    {
        gmx_exit_on_fatal_error(ExitType_Abort, 1);
    }
#endif

    return (wanted == compiled);
}

/*! \endcond */

} // namespace gmx
