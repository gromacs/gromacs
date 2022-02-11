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
#ifndef GMX_SIMD_SUPPORT_H
#define GMX_SIMD_SUPPORT_H


/*! \libinternal \file
 *
 * \brief Functions to query compiled and supported SIMD architectures
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \inlibraryapi
 * \ingroup module_simd
 */

#include "gromacs/hardware/cpuinfo.h"

namespace gmx
{

/*! \cond libapi */

/*! \brief Enumerated options for SIMD architectures */
enum class SimdType
{
    None,          //!< Disable all SIMD support
    Reference,     //!< Gromacs reference software SIMD
    Generic,       //!< Placeholder for future support for gcc generic SIMD
    X86_Sse2,      //!< SSE2
    X86_Sse4_1,    //!< SSE4.1
    X86_Avx128Fma, //!< 128-bit Avx with FMA (Amd)
    X86_Avx,       //!< 256-bit Avx
    X86_Avx2,      //!< AVX2
    X86_Avx2_128,  //!< 128-bit AVX2, better than 256-bit for AMD Ryzen
    X86_Avx512,    //!< AVX_512
    X86_Avx512Knl, //!< AVX_512_KNL
    Arm_NeonAsimd, //!< 64-bit ARM AArch64 Advanced SIMD
    Arm_Sve,       //!< ARM Scalable Vector Extensions
    Ibm_Vsx        //!< IBM VSX SIMD (Power7 and later)
};

/*! \libinternal \brief Return a string with the name of a SIMD type
 *
 *  \param s  SIMD type to turn into string
 */
const std::string& simdString(SimdType s);

/*! \libinternal \brief Return the SIMD type that would fit this hardware best */
SimdType simdSuggested(const CpuInfo& c);

/*! \libinternal \brief Return the SIMD type the library was compiled with */
SimdType simdCompiled();

/*! \libinternal \brief Check if binary was compiled with the provided SIMD type
 *
 *  \param s              SIMD type to query. If this matches the suggested type
 *                        for this cpu, the routine returns quietly.
 *  \param log            If not nullptr, statistics will be printed to the file.
 *                        If we do not have a match there will also be a warning.
 *  \param warnToStdErr   If true, warnings will also be printed to stderr.
 */
bool simdCheck(SimdType s, FILE* log, bool warnToStdErr);

/*! \endcond */

} // namespace gmx


#endif // GMX_SIMD_SUPPORT_H
