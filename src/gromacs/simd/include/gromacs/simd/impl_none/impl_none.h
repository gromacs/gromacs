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

#ifndef GMX_SIMD_IMPL_NONE_H
#define GMX_SIMD_IMPL_NONE_H

/* No SIMD implementation - assign 0 to all defines */
#define GMX_SIMD 0
#define GMX_SIMD_HAVE_FLOAT 0
#define GMX_SIMD_HAVE_DOUBLE 0
#define GMX_SIMD_HAVE_LOADU 0
#define GMX_SIMD_HAVE_STOREU 0
#define GMX_SIMD_HAVE_LOGICAL 0
#define GMX_SIMD_HAVE_FMA 0
#define GMX_SIMD_HAVE_FINT32_EXTRACT 0
#define GMX_SIMD_HAVE_FINT32_LOGICAL 0
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS 0
#define GMX_SIMD_HAVE_DINT32_EXTRACT 0
#define GMX_SIMD_HAVE_DINT32_LOGICAL 0
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS 0
#define GMX_SIMD4_HAVE_FLOAT 0
#define GMX_SIMD4_HAVE_DOUBLE 0
#define GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT 0
#define GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE 0

#undef GMX_SIMD_FLOAT_WIDTH
#undef GMX_SIMD_DOUBLE_WIDTH
#undef GMX_SIMD_FINT32_WIDTH
#undef GMX_SIMD_DINT32_WIDTH
#undef GMX_SIMD4_WIDTH
#define GMX_SIMD_ALIGNMENT 8 // 1*double
#undef GMX_SIMD_RSQRT_BITS
#undef GMX_SIMD_RCP_BITS

#endif /* GMX_SIMD_IMPL_NONE_H */
