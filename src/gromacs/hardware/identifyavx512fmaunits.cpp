/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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

#ifdef GMX_IDENTIFY_AVX512_FMA_UNITS_STANDALONE
#define GMX_ENABLE_AVX512_TESTS 1
#else
#include "config.h"
#endif

#if GMX_ENABLE_AVX512_TESTS
#include <immintrin.h>

#ifdef _MSC_VER
#include <intrin.h>
#endif
#endif // GMX_ENABLE_AVX512_TESTS

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

#if GMX_ENABLE_AVX512_TESTS
// Use a local routine to read the timestep counter just on x86 to avoid dependence
// on the Gromacs cycle counter module.
uint64_t rdtscp(void)
{
#ifdef MSC_VER
    unsigned int ui;
    return static_cast<uint64_t>(__rdtscp(&ui));
#else
    uint32_t low;
    uint32_t high;

    __asm__ __volatile__("rdtscp" : "=a" (low), "=d" (high) :: "ecx" );
    return (static_cast<uint64_t>(high) << 32) | low;
#endif
}

/*\ brief Loop over mixed FMA and shuffle AVX512 instructions
 *
 * This function executes a meaningless loop that includes both
 * FMA and shuffle instructions from the AVX512 instruction set.
 * We need a bit of complex logic to make sure it cannot be
 * optimized away by the compiler.
 *
 * \param loopCount  Number of iterations. Each iteration will
 *                   execute 12 FMA and 12 shuffle instructions.
 * \param seed       A double-precision number between 0 and 1.
 *                   To be really certain the loop is not optimized
 *                   away, you should use some timing-related
 *                   function to create this seed at runtime.
 * \return Meaningless floating-point number. Make sure you
 *         add this number to some variable and conditionally
 *         issue a print statement e.g. if it is negative
 *         (which should not happen), again to make sure the loop
 *         cannot be optimized away.
 */
double
executeFmaAndShuffleLoop(int     loopCount,
                         double  seed)
{
    // Make sure all variables are different to avoid gcc optimizing them away
    __m512d   d0   = _mm512_set1_pd(1.0-0.01*seed);
    __m512d   d1   = _mm512_set1_pd(1.0-0.02*seed);
    __m512d   d2   = _mm512_set1_pd(1.0-0.03*seed);
    __m512d   d3   = _mm512_set1_pd(1.0-0.04*seed);
    __m512d   d4   = _mm512_set1_pd(1.0-0.05*seed);
    __m512d   d5   = _mm512_set1_pd(1.0-0.06*seed);
    __m512d   d6   = _mm512_set1_pd(1.0-0.07*seed);
    __m512d   d7   = _mm512_set1_pd(1.0-0.08*seed);
    __m512d   d8   = _mm512_set1_pd(1.0-0.09*seed);
    __m512d   d9   = _mm512_set1_pd(1.0-0.10*seed);
    __m512d   d10  = _mm512_set1_pd(1.0-0.11*seed);
    __m512d   d11  = _mm512_set1_pd(1.0-0.12*seed);
    __m512d   eps  = _mm512_set1_pd(1e-6);
    __m512i   i0   = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    __m512i   i1   = _mm512_set_epi32(0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1);
    __m512i   i2   = _mm512_set_epi32(1, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2);
    __m512i   i3   = _mm512_set_epi32(2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3);
    __m512i   i4   = _mm512_set_epi32(3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4);
    __m512i   i5   = _mm512_set_epi32(4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5);
    __m512i   i6   = _mm512_set_epi32(5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6);
    __m512i   i7   = _mm512_set_epi32(7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8);
    __m512i   i8   = _mm512_set_epi32(8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9);
    __m512i   i9   = _mm512_set_epi32(9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10);
    __m512i   i10  = _mm512_set_epi32(10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11);
    __m512i   i11  = _mm512_set_epi32(11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12);
    __m512i   idx  = _mm512_set_epi32(12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13);
    __mmask16 mask = static_cast<uint16_t>(0xffff);

    for (int i = 0; i < loopCount; i++)
    {
        d0  = _mm512_fmadd_pd(d0, d0, eps);
        d1  = _mm512_fmadd_pd(d1, d1, eps);
        d2  = _mm512_fmadd_pd(d2, d2, eps);
        d3  = _mm512_fmadd_pd(d3, d3, eps);
        d4  = _mm512_fmadd_pd(d4, d4, eps);
        d5  = _mm512_fmadd_pd(d5, d5, eps);
        d6  = _mm512_fmadd_pd(d6, d6, eps);
        d7  = _mm512_fmadd_pd(d7, d7, eps);
        d8  = _mm512_fmadd_pd(d8, d8, eps);
        d9  = _mm512_fmadd_pd(d9, d9, eps);
        d10 = _mm512_fmadd_pd(d10, d10, eps);
        d11 = _mm512_fmadd_pd(d11, d11, eps);
        // plain permutevar is not yet available in gcc-6.4
        i0  = _mm512_maskz_permutexvar_epi32(mask, idx, i0);
        i1  = _mm512_maskz_permutexvar_epi32(mask, idx, i1);
        i2  = _mm512_maskz_permutexvar_epi32(mask, idx, i2);
        i3  = _mm512_maskz_permutexvar_epi32(mask, idx, i3);
        i4  = _mm512_maskz_permutexvar_epi32(mask, idx, i4);
        i5  = _mm512_maskz_permutexvar_epi32(mask, idx, i5);
        i6  = _mm512_maskz_permutexvar_epi32(mask, idx, i6);
        i7  = _mm512_maskz_permutexvar_epi32(mask, idx, i7);
        i8  = _mm512_maskz_permutexvar_epi32(mask, idx, i8);
        i9  = _mm512_maskz_permutexvar_epi32(mask, idx, i9);
        i10 = _mm512_maskz_permutexvar_epi32(mask, idx, i10);
        i11 = _mm512_maskz_permutexvar_epi32(mask, idx, i11);
    }

    // Make sure we use all variables in the loop to return a result
    i0  = _mm512_add_epi32(i0, i1);
    i2  = _mm512_add_epi32(i2, i3);
    i4  = _mm512_add_epi32(i4, i5);
    i6  = _mm512_add_epi32(i6, i7);
    i8  = _mm512_add_epi32(i8, i9);
    i10 = _mm512_add_epi32(i10, i11);
    i0  = _mm512_add_epi32(i0, i2);
    i4  = _mm512_add_epi32(i4, i6);
    i8  = _mm512_add_epi32(i8, i10);
    i0  = _mm512_add_epi32(i0, i4);
    i0  = _mm512_add_epi32(i0, i8);

    d0 = _mm512_fmadd_pd(d0, d1, d2);
    d3 = _mm512_fmadd_pd(d3, d4, d5);
    d6 = _mm512_fmadd_pd(d6, d7, d8);
    d9 = _mm512_fmadd_pd(d9, d10, d11);
    d0 = _mm512_add_pd(d0, d3);
    d6 = _mm512_add_pd(d6, d9);
    d0 = _mm512_add_pd(d0, d6);

    double data[8];
    int    idata[16];
    _mm512_storeu_pd(data, d0);
    _mm512_storeu_si512(idata, i0);

    double d = 0;

    for (int i = 0; i < 8; i++)
    {
        d += data[i] * idata[2*i] * idata[2*i+1];
    }

    return d;
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
 * \param seed       A double-precision number between 0 and 1.
 *                   To be really certain the loop is not optimized
 *                   away, you should use some timing-related
 *                   function to create this seed at runtime.
 * \return Meaningless floating-point number. Make sure you
 *         add this number to some variable and conditionally
 *         issue a print statement e.g. if it is negative
 *         (which should not happen), again to make sure the loop
 *         cannot be optimized away.
 */
double
executeFmaOnlyLoop(int    loopCount,
                   double seed)
{
    // Make sure all variables are different to avoid gcc optimizing them away
    __m512d d0  = _mm512_set1_pd(1.0-0.01*seed);
    __m512d d1  = _mm512_set1_pd(1.0-0.02*seed);
    __m512d d2  = _mm512_set1_pd(1.0-0.03*seed);
    __m512d d3  = _mm512_set1_pd(1.0-0.04*seed);
    __m512d d4  = _mm512_set1_pd(1.0-0.05*seed);
    __m512d d5  = _mm512_set1_pd(1.0-0.06*seed);
    __m512d d6  = _mm512_set1_pd(1.0-0.07*seed);
    __m512d d7  = _mm512_set1_pd(1.0-0.08*seed);
    __m512d d8  = _mm512_set1_pd(1.0-0.09*seed);
    __m512d d9  = _mm512_set1_pd(1.0-0.10*seed);
    __m512d d10 = _mm512_set1_pd(1.0-0.11*seed);
    __m512d d11 = _mm512_set1_pd(1.0-0.12*seed);
    __m512d eps = _mm512_set1_pd(1e-6);

    for (int i = 0; i < loopCount; i++)
    {
        d0  = _mm512_fmadd_pd(d0, d0, eps);
        d1  = _mm512_fmadd_pd(d1, d1, eps);
        d2  = _mm512_fmadd_pd(d2, d2, eps);
        d3  = _mm512_fmadd_pd(d3, d3, eps);
        d4  = _mm512_fmadd_pd(d4, d4, eps);
        d5  = _mm512_fmadd_pd(d5, d5, eps);
        d6  = _mm512_fmadd_pd(d6, d6, eps);
        d7  = _mm512_fmadd_pd(d7, d7, eps);
        d8  = _mm512_fmadd_pd(d8, d8, eps);
        d9  = _mm512_fmadd_pd(d9, d9, eps);
        d10 = _mm512_fmadd_pd(d10, d10, eps);
        d11 = _mm512_fmadd_pd(d11, d11, eps);
    }

    // Make sure we use all variables in the loop to return a result
    d0 = _mm512_fmadd_pd(d0, d1, d2);
    d3 = _mm512_fmadd_pd(d3, d4, d5);
    d6 = _mm512_fmadd_pd(d6, d7, d8);
    d9 = _mm512_fmadd_pd(d9, d10, d11);
    d0 = _mm512_add_pd(d0, d3);
    d6 = _mm512_add_pd(d6, d9);
    d0 = _mm512_add_pd(d0, d6);

    double data[8];

    _mm512_storeu_pd(data, d0);

    double d = 0;

    for (int i = 0; i < 8; i++)
    {
        d += data[i];
    }
    return d;
}

int
checkDualAvx512FmaUnits()
{
    uint64_t timeFmaAndShuf = 1e9;             // Large value
    uint64_t timeFmaOnly    = 1e9;             // Large value
    double   dummy;
    double   seed = (rdtscp() & 0xff) / 256.0; // Create an unpredictable small number between 0 and 1

    // Make sure the CPU is in AVX512 mode by executing a fairly long loop
    dummy = executeFmaOnlyLoop(100000, seed);

    // Execute the mixed FMA/shuffle loop three times
    for (int i = 0; i < 3; i++)
    {
        uint64_t start = rdtscp();
        dummy += executeFmaAndShuffleLoop(1000, seed);
        uint64_t res = rdtscp() - start;
        timeFmaAndShuf = std::min(timeFmaAndShuf, res);
    }

    // Execute the FMA-only loop three times
    for (int i = 0; i < 3; i++)
    {
        uint64_t start = rdtscp();
        dummy += executeFmaOnlyLoop(1000, seed);
        uint64_t res = rdtscp() - start;
        timeFmaOnly = std::min(timeFmaOnly, res);
    }

    // Dummy can never be negative, but by using it in the
    // conditional it cannot be optimized away.
    return (timeFmaAndShuf > 1.5 * timeFmaOnly || dummy < 0);
}


#endif  // GMX_ENABLE_AVX512_TESTS

std::mutex initMutex;

}      // namespace anonymous

int
identifyAvx512FmaUnits()
{
    static bool initialized = false;
    static int  result      = false;

    if (!initialized)
    {
        initMutex.lock();
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
#if GMX_ENABLE_AVX512_TESTS
                result = checkDualAvx512FmaUnits() ? 2 : 1;
#else
                result = -1; // Cannot run the tests
#endif
            }
            else
            {
                result = 0; // Not AVX-512 hardware
            }

        }
        initMutex.unlock();
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
