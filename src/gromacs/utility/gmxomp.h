/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares OpenMP wrappers to avoid conditional compilation.
 *
 * This module defines wrappers for OpenMP API functions and enables compiling
 * code without conditional compilation even when OpenMP is turned off in the
 * build system.
 * Therefore, OpenMP API functions should always be used through these wrappers
 * and omp.h should never be directly included.  Instead, this header should be
 * used whenever OpenMP API functions are needed.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_OMP_H
#define GMX_UTILITY_OMP_H

#include "config.h"

#include <stdio.h>

#ifndef GMX_NATIVE_WINDOWS
/* Ugly hack because the openmp implementation below hacks into the SIMD
 * settings to decide when to use _mm_pause(). This should eventually be
 * changed into proper detection of the intrinsics uses, not SIMD.
 */
#if (defined GMX_SIMD_X86_SSE2) || (defined GMX_SIMD_X86_SSE4_1) || \
    (defined GMX_SIMD_X86_AVX_128_FMA) || (defined GMX_SIMD_X86_AVX_256) || \
    (defined GMX_SIMD_X86_AVX2_256)
#    include <xmmintrin.h>
#endif
#else
#include <windows.h>
#endif

#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*! \addtogroup module_utility
 * \{
 */

/*! \brief
 * Returns an integer equal to or greater than the number of threads
 * that would be available if a parallel region without num_threads were
 * defined at that point in the code.
 *
 * Acts as a wrapper for omp_get_max_threads().
 */
int gmx_omp_get_max_threads(void);

/*! \brief
 * Returns the number of processors available when the function is called.
 *
 * Acts as a wrapper around omp_get_num_procs().
 */
int gmx_omp_get_num_procs(void);

/*! \brief
 * Returns the thread number of the thread executing within its thread team.
 *
 * Acts as a wrapper for omp_get_thread_num().
 */
int gmx_omp_get_thread_num(void);

/*! \brief
 * Sets the number of threads in subsequent parallel regions, unless overridden
 * by a num_threads clause.
 *
 * Acts as a wrapper for omp_set_num_threads().
 */
void gmx_omp_set_num_threads(int num_threads);

/*! \brief
 * Check for externally set thread affinity to avoid conflicts with \Gromacs
 * internal setting.
 *
 * \param[out] message  Receives the message to be shown to the user.
 * \returns `true` if we can set thread affinity ourselves.
 *
 * While GNU OpenMP does not set affinity by default, the Intel OpenMP library
 * does.  This conflicts with the internal affinity (especially thread-MPI)
 * setting, results in incorrectly locked threads, and causes dreadful performance.
 *
 * The KMP_AFFINITY environment variable is used by Intel, GOMP_CPU_AFFINITY
 * by the GNU compilers (Intel also honors it well).  If any of the variables
 * is set, we should honor it and disable the internal pinning.
 * When using Intel OpenMP, we will disable affinity if the user did not set it
 * manually through one of the aforementioned environment variables.
 *
 * Note that the Intel OpenMP affinity disabling will only take effect if this
 * function is called before the OpenMP library gets initialized, which happens
 * when the first call is made into a compilation unit that contains OpenMP
 * pragmas.
 *
 * If this function returns `false`, the caller is responsible to disable the
 * pinning, show the message from \p *message to the user, and free the memory
 * allocated for \p *message.
 * If the return value is `true`, \p *message is NULL.
 */
gmx_bool gmx_omp_check_thread_affinity(char **message);

/*! \brief
 * Pause for use in a spin-wait loop.
 */
static gmx_inline void gmx_pause()
{
#ifndef _MSC_VER
    /* Ugly hack because the openmp implementation below hacks into the SIMD
     * settings to decide when to use _mm_pause(). This should eventually be
     * changed into proper detection of the intrinsics uses, not SIMD.
     */
#if ((defined GMX_SIMD_X86_SSE2) || (defined GMX_SIMD_X86_SSE4_1) || \
    (defined GMX_SIMD_X86_AVX_128_FMA) || (defined GMX_SIMD_X86_AVX_256) || \
    (defined GMX_SIMD_X86_AVX2_256)) && !defined(__MINGW32__)
    /* Replace with tbb::internal::atomic_backoff when/if we use TBB */
    _mm_pause();
#elif defined __MIC__
    _mm_delay_32(32);
#else
    /* No wait for unknown architecture */
#endif
#else
    YieldProcessor();
#endif
}

/*! \} */

#ifdef __cplusplus
}
#endif

#endif
