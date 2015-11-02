/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_SPARC64_HPC_ACE_COMMON_H
#define GMX_SIMD_IMPL_SPARC64_HPC_ACE_COMMON_H

/* Capability definitions for Sparc64 HPC-ACE */
/* HPC-ACE is actually double-only on the register level, but we also implement
 * a single-precision interface where we only offer single-precision accuracy
 * in math functions - this can save quite a few cycles.
 */
#define GMX_SIMD                                 1
#define GMX_SIMD_HAVE_FLOAT                      1
#define GMX_SIMD_HAVE_DOUBLE                     1
#define GMX_SIMD_HAVE_LOADU                      0
#define GMX_SIMD_HAVE_STOREU                     0
#define GMX_SIMD_HAVE_LOGICAL                    1
#define GMX_SIMD_HAVE_FMA                        1
#define GMX_SIMD_HAVE_FRACTION                   0
#define GMX_SIMD_HAVE_FINT32                     1
#define GMX_SIMD_HAVE_FINT32_EXTRACT             1
#define GMX_SIMD_HAVE_FINT32_LOGICAL             1
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS         0
#define GMX_SIMD_HAVE_DINT32                     1
#define GMX_SIMD_HAVE_DINT32_EXTRACT             1
#define GMX_SIMD_HAVE_DINT32_LOGICAL             1
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS         0
#define GMX_SIMD4_HAVE_FLOAT                     0
#define GMX_SIMD4_HAVE_DOUBLE                    0

/* Implementation details */
#define GMX_SIMD_FLOAT_WIDTH                     2
#define GMX_SIMD_DOUBLE_WIDTH                    2
#define GMX_SIMD_FINT32_WIDTH                    2
#define GMX_SIMD_DINT32_WIDTH                    2
/* No SIMD4 support in either single or double */
#define GMX_SIMD_RSQRT_BITS                     10
#define GMX_SIMD_RCP_BITS                        9

#endif /* GMX_SIMD_IMPL_SPARC64_HPC_ACE_COMMON_H */
