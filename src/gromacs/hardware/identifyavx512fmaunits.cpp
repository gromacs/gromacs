/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief Implements a routine to check the number of AVX512 fma units
 *
 * Just as the CpuInfo code, we need to be able to compile this file in stand-alone mode
 * to set the SIMD acceleration and similar things during CMake configuration.
 */

#ifndef GMX_IDENTIFY_AVX512_FMA_UNITS_STANDALONE
#include "gmxpre.h"
#endif

#include "identifyavx512fmaunits.h"

#ifndef GMX_IDENTIFY_AVX512_FMA_UNITS_STANDALONE
#include "config.h"
#endif

#include <cstdint>
#include <cstdio>

#include <algorithm>
#include <mutex>

#ifndef GMX_IDENTIFY_AVX512_FMA_UNITS_STANDALONE
#include "gromacs/hardware/cpuinfo.h"
#endif

namespace gmx
{

namespace
{

#if GMX_X86_GCC_INLINE_ASM && SIMD_AVX_512_CXX_SUPPORTED
/*\ brief Loop over mixed FMA and shuffle AVX512 instructions
 *
 * This function executes a meaningless loop that includes both
 * FMA and shuffle instructions from the AVX512 instruction set.
 * We need a bit of complex logic to make sure it cannot be
 * optimized away by the compiler.
 *
 * \param loopCount  Number of iterations. Each iteration will
 *                   execute 12 FMA and 12 shuffle instructions.
 * \return Number of cycles used for the loop.
 */
uint64_t
timeFmaAndShuffleLoop(uint64_t loopCount)
{
    uint64_t cycles;
    // Unfortunately we need to resort to inline ASM since we are
    // making a choice based on timing, and without efficient optimization
    // (e.g. when doing debugging) the usual intrinsics are often implemented
    // as independent load/store operations, which completely screws up timing.
    __asm__ __volatile__("\tvpxord %%zmm0, %%zmm0, %%zmm0\n"
                         "\tvmovaps %%zmm0, %%zmm1\n"
                         "\tvmovaps %%zmm0, %%zmm2\n"
                         "\tvmovaps %%zmm0, %%zmm3\n"
                         "\tvmovaps %%zmm0, %%zmm4\n"
                         "\tvmovaps %%zmm0, %%zmm5\n"
                         "\tvmovaps %%zmm0, %%zmm6\n"
                         "\tvmovaps %%zmm0, %%zmm7\n"
                         "\tvmovaps %%zmm0, %%zmm8\n"
                         "\tvmovaps %%zmm0, %%zmm9\n"
                         "\tvmovaps %%zmm0, %%zmm10\n"
                         "\tvmovaps %%zmm0, %%zmm11\n"
                         "\tvpxord %%zmm12, %%zmm12, %%zmm12\n"
                         "\tvmovaps %%zmm12, %%zmm13\n"
                         "\tvmovaps %%zmm12, %%zmm14\n"
                         "\tvmovaps %%zmm12, %%zmm15\n"
                         "\tvmovaps %%zmm12, %%zmm16\n"
                         "\tvmovaps %%zmm12, %%zmm17\n"
                         "\tvmovaps %%zmm12, %%zmm18\n"
                         "\tvmovaps %%zmm12, %%zmm19\n"
                         "\tvmovaps %%zmm12, %%zmm20\n"
                         "\tvmovaps %%zmm12, %%zmm21\n"
                         "\tvmovaps %%zmm12, %%zmm22\n"
                         "\tvmovaps %%zmm12, %%zmm23\n"
                         "\tvmovaps %%zmm12, %%zmm30\n"
                         "\trdtscp\n"
                         "\tsalq $32, %%rdx\n"
                         "\tmovl %%eax, %%eax\n"
                         "\tmovq %%rdx, %%rbx\n"
                         "\torq %%rax, %%rbx\n"
                         "\tmovq %1, %%rdx\n"
                         "1:\n"
                         "\tvfmadd231pd %%zmm0, %%zmm0, %%zmm0\n"
                         "\tvfmadd231pd %%zmm1, %%zmm1, %%zmm1\n"
                         "\tvfmadd231pd %%zmm2, %%zmm2, %%zmm2\n"
                         "\tvfmadd231pd %%zmm3, %%zmm3, %%zmm3\n"
                         "\tvfmadd231pd %%zmm4, %%zmm4, %%zmm4\n"
                         "\tvfmadd231pd %%zmm5, %%zmm5, %%zmm5\n"
                         "\tvfmadd231pd %%zmm6, %%zmm6, %%zmm6\n"
                         "\tvfmadd231pd %%zmm7, %%zmm7, %%zmm7\n"
                         "\tvfmadd231pd %%zmm8, %%zmm8, %%zmm8\n"
                         "\tvfmadd231pd %%zmm9, %%zmm9, %%zmm9\n"
                         "\tvfmadd231pd %%zmm10, %%zmm10, %%zmm10\n"
                         "\tvfmadd231pd %%zmm11, %%zmm11, %%zmm11\n"
                         "\tvpermd %%zmm30, %%zmm30, %%zmm12\n"
                         "\tvpermd %%zmm30, %%zmm30, %%zmm13\n"
                         "\tvpermd %%zmm30, %%zmm30, %%zmm14\n"
                         "\tvpermd %%zmm30, %%zmm30, %%zmm15\n"
                         "\tvpermd %%zmm30, %%zmm30, %%zmm16\n"
                         "\tvpermd %%zmm30, %%zmm30, %%zmm17\n"
                         "\tvpermd %%zmm30, %%zmm30, %%zmm18\n"
                         "\tvpermd %%zmm30, %%zmm30, %%zmm19\n"
                         "\tvpermd %%zmm30, %%zmm30, %%zmm20\n"
                         "\tvpermd %%zmm30, %%zmm30, %%zmm21\n"
                         "\tvpermd %%zmm30, %%zmm30, %%zmm22\n"
                         "\tvpermd %%zmm30, %%zmm30, %%zmm23\n"
                         "\tdec %%rdx\n"
                         "\tjg 1b\n"
                         "\trdtscp\n"
                         "\tsalq $32, %%rdx\n"
                         "\tmovl %%eax, %%eax\n"
                         "\torq %%rax, %%rdx\n"
                         "\tsubq %%rbx, %%rdx\n"
                         "\tmovq %%rdx, %0\n"
                         : "=r" (cycles) : "r" (loopCount)
                         : "rax", "rbx", "rcx", "rdx", "zmm0", "zmm1", "zmm2", "zmm3",
                         "zmm4", "zmm5", "zmm6", "zmm7", "zmm8", "zmm9", "zmm10",
                         "zmm11", "zmm12", "zmm13", "zmm14", "zmm15", "zmm16", "zmm17",
                         "zmm18", "zmm19", "zmm20", "zmm21", "zmm22", "zmm23", "zmm30");

    return cycles;
}

/*\ brief Loop over FMA AVX512 instructions
 *
 * This function executes a meaningless loop that includes only
 * FMA instructions from the AVX512 instruction set.
 * We need a bit of complex logic to make sure it cannot be
 * optimized away by the compiler.
 *
 * \param loopCount  Number of iterations. Each iteration will
 *                   execute 12 FMA instructions.
 * \return Number of cycles used for the loop.
 */
uint64_t
timeFmaOnlyLoop(uint64_t loopCount)
{
    uint64_t cycles;
    // Unfortunately we need to resort to inline ASM since we are
    // making a choice based on timing, and without efficient optimization
    // (e.g. when doing debugging) the usual intrinsics are often implemented
    // as independent load/store operations, which completely screws up timing.
    __asm__ __volatile__("\tvpxord %%zmm0, %%zmm0, %%zmm0\n"
                         "\tvmovaps %%zmm0, %%zmm1\n"
                         "\tvmovaps %%zmm0, %%zmm2\n"
                         "\tvmovaps %%zmm0, %%zmm3\n"
                         "\tvmovaps %%zmm0, %%zmm4\n"
                         "\tvmovaps %%zmm0, %%zmm5\n"
                         "\tvmovaps %%zmm0, %%zmm6\n"
                         "\tvmovaps %%zmm0, %%zmm7\n"
                         "\tvmovaps %%zmm0, %%zmm8\n"
                         "\tvmovaps %%zmm0, %%zmm9\n"
                         "\tvmovaps %%zmm0, %%zmm10\n"
                         "\tvmovaps %%zmm0, %%zmm11\n"
                         "\trdtscp\n"
                         "\tsalq $32, %%rdx\n"
                         "\tmovl %%eax, %%eax\n"
                         "\tmovq %%rdx, %%rbx\n"
                         "\torq %%rax, %%rbx\n"
                         "\tmovq %1, %%rdx\n"
                         "1:\n"
                         "\tvfmadd231pd %%zmm0, %%zmm0, %%zmm0\n"
                         "\tvfmadd231pd %%zmm1, %%zmm1, %%zmm1\n"
                         "\tvfmadd231pd %%zmm2, %%zmm2, %%zmm2\n"
                         "\tvfmadd231pd %%zmm3, %%zmm3, %%zmm3\n"
                         "\tvfmadd231pd %%zmm4, %%zmm4, %%zmm4\n"
                         "\tvfmadd231pd %%zmm5, %%zmm5, %%zmm5\n"
                         "\tvfmadd231pd %%zmm6, %%zmm6, %%zmm6\n"
                         "\tvfmadd231pd %%zmm7, %%zmm7, %%zmm7\n"
                         "\tvfmadd231pd %%zmm8, %%zmm8, %%zmm8\n"
                         "\tvfmadd231pd %%zmm9, %%zmm9, %%zmm9\n"
                         "\tvfmadd231pd %%zmm10, %%zmm10, %%zmm10\n"
                         "\tvfmadd231pd %%zmm11, %%zmm11, %%zmm11\n"
                         "\tdec %%rdx\n"
                         "\tjg 1b\n"
                         "\trdtscp\n"
                         "\tsalq $32, %%rdx\n"
                         "\tmovl %%eax, %%eax\n"
                         "\torq %%rax, %%rdx\n"
                         "\tsubq %%rbx, %%rdx\n"
                         "\tmovq %%rdx, %0\n"
                         : "=r" (cycles) : "r" (loopCount)
                         : "rax", "rbx", "rcx", "rdx", "zmm0", "zmm1", "zmm2", "zmm3",
                         "zmm4", "zmm5", "zmm6", "zmm7", "zmm8", "zmm9", "zmm10", "zmm11");

    return cycles;
}

int
checkDualAvx512FmaUnits()
{
    uint64_t timeFmaAndShuf = 1e9;             // Large value

    // Make sure the CPU is in AVX512 mode by executing a fairly long loop.
    // Use the return value to make sure it is not optimized away. Later invocations
    // use fewer iterations, so they should always be faster.
    uint64_t timeFmaOnly    = timeFmaOnlyLoop(100000);

    // Execute the loops three times
    for (int i = 0; i < 3; i++)
    {
        timeFmaAndShuf = std::min(timeFmaAndShuf, timeFmaAndShuffleLoop(1000) );
        timeFmaOnly    = std::min(timeFmaOnly, timeFmaOnlyLoop(1000) );
    }

    return (timeFmaAndShuf > 1.5 * timeFmaOnly);
}

#endif  // GMX_X86_GCC_INLINE_ASM && SIMD_AVX_512_CXX_SUPPORTED


/*! \brief Mutex to guard the execution of the timing test
 *
 * We only execute the test once, and return the saved result
 * on subsequent calls.
 */
std::mutex initMutex;

}      // namespace anonymous

int
identifyAvx512FmaUnits()
{
    static bool initialized = false;
    static int  result      = false;

    if (!initialized)
    {
        std::lock_guard<std::mutex>  lock(initMutex);

        if (!initialized)
        {
            // For the standalone test binary we assume it will
            // only be executed on AVX512 hardware, but for the
            // library version we check the hardware support.
#ifdef GMX_IDENTIFY_AVX512_FMA_UNITS_STANDALONE
            bool haveAvx512Hardware = true;
#else
            bool haveAvx512Hardware = CpuInfo::detect().feature(CpuInfo::Feature::X86_Avx512F);
#endif

            if (haveAvx512Hardware)
            {
#if GMX_X86_GCC_INLINE_ASM && SIMD_AVX_512_CXX_SUPPORTED
                result = checkDualAvx512FmaUnits() ? 2 : 1;
#else
                result = -1; // Cannot run the tests
#endif
            }
            else
            {
                result = 0; // Not AVX-512 hardware
            }
            initialized = true;
        }
    }
    return result;
}

} // namespace gmx

#ifdef GMX_IDENTIFY_AVX512_FMA_UNITS_STANDALONE
int
main()
{
    printf("%d\n", gmx::identifyAvx512FmaUnits());
    return 0;
}
#endif
