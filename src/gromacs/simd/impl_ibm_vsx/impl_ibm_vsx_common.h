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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VSX_COMMON_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VSX_COMMON_H

/* Since we've had to use quite a few compiler and endian defines in this code
 * it might not 'just work' when e.g. clang provides VSX support. If you are
 * reading this because you saw the error below for a new compiler, try removing
 * the checks, but then make sure you run 'make test'.
 */
#if !(defined __GNUC__) && !(defined __ibmxl__) && !(defined __xlC__)
#    error VSX acceleration is very compiler-dependent, and only tested for gcc & xlc.
#endif

/* Capability definitions for IBM VSX */
#define GMX_SIMD                               1
#define GMX_SIMD_HAVE_FLOAT                    1
#define GMX_SIMD_HAVE_LOADU                    1
#define GMX_SIMD_HAVE_STOREU                   1
#define GMX_SIMD_HAVE_LOGICAL                  1
#define GMX_SIMD_HAVE_FMA                      1
#define GMX_SIMD_HAVE_FRACTION                 0
#define GMX_SIMD_HAVE_FINT32                   1
#define GMX_SIMD_HAVE_FINT32_EXTRACT           1
#define GMX_SIMD_HAVE_FINT32_LOGICAL           1
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS       1

/* With GCC, only version 4.9 or later supports all parts of double precision VSX.
 * We check explicitly for xlc, since that compiler appears to like pretending it is gcc,
 * but there double precision seems to work fine.
 */
#if defined(__ibmxl__) || defined(__xlC__) || !(defined(__GNUC__) && ((__GNUC__ < 4) || ((__GNUC__ == 4) && (__GNUC_MINOR < 9))))
#    define GMX_SIMD_HAVE_DOUBLE               1
#    define GMX_SIMD_HAVE_DINT32               1
#    define GMX_SIMD_HAVE_DINT32_EXTRACT       1
#    define GMX_SIMD_HAVE_DINT32_LOGICAL       1
#    define GMX_SIMD_HAVE_DINT32_ARITHMETICS   1
#else
#    define GMX_SIMD_HAVE_DOUBLE               0
#    define GMX_SIMD_HAVE_DINT32               0
#    define GMX_SIMD_HAVE_DINT32_EXTRACT       0
#    define GMX_SIMD_HAVE_DINT32_LOGICAL       0
#    define GMX_SIMD_HAVE_DINT32_ARITHMETICS   0
#endif

#define GMX_SIMD4_HAVE_FLOAT                   1
#define GMX_SIMD4_HAVE_DOUBLE                  0

/* Implementation details */
#define GMX_SIMD_FLOAT_WIDTH                   4
#define GMX_SIMD_DOUBLE_WIDTH                  2
#define GMX_SIMD_FINT32_WIDTH                  4
#define GMX_SIMD_DINT32_WIDTH                  2
#define GMX_SIMD4_WIDTH                        4
#define GMX_SIMD_RSQRT_BITS                   14
#define GMX_SIMD_RCP_BITS                     14

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_VSX_COMMON_H */
