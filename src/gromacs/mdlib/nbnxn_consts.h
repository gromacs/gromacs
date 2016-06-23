/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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

#ifndef _nbnxn_consts_h
#define _nbnxn_consts_h

#ifdef __cplusplus
extern "C" {
#endif


/* With CPU kernels the i-cluster size is always 4 atoms.
 * With x86 SIMD the j-cluster size can be 2, 4 or 8, otherwise 4.
 */
#define NBNXN_CPU_CLUSTER_I_SIZE       4

#define NBNXN_CPU_CLUSTER_I_SIZE_2LOG  2

// Lower limit for square interaction distances in nonbonded kernels.
// For smaller values we will overflow when calculating r^-1 or r^-12, but
// to keep it simple we always apply the limit from the tougher r^-12 condition.
#if GMX_DOUBLE
// Some double precision SIMD architectures use single precision in the first
// step, so although the double precision criterion would allow smaller rsq,
// we need to stay in single precision with some margin for the N-R iterations.
#define NBNXN_MIN_RSQ         1.0e-36
#else
// The worst intermediate value we might evaluate is r^-12, which
// means we should ensure r^2 stays above pow(GMX_FLOAT_MAX,-1.0/6.0)*1.01 (some margin)
#define NBNXN_MIN_RSQ         3.82e-07f  // r > 6.2e-4
#endif


/* Cluster-pair Interaction masks for 4xN and 2xNN kernels.
 * Bit i*CJ_SIZE + j tells if atom i and j interact.
 */
/* All interaction mask is the same for all kernels */
#define NBNXN_INTERACTION_MASK_ALL        0xffffffffU
/* 4x4 kernel diagonal mask */
#define NBNXN_INTERACTION_MASK_DIAG       0x08ceU
/* 4x2 kernel diagonal masks */
#define NBNXN_INTERACTION_MASK_DIAG_J2_0  0x0002U
#define NBNXN_INTERACTION_MASK_DIAG_J2_1  0x002fU
/* 4x8 kernel diagonal masks */
#define NBNXN_INTERACTION_MASK_DIAG_J8_0  0xf0f8fcfeU
#define NBNXN_INTERACTION_MASK_DIAG_J8_1  0x0080c0e0U


#ifdef __cplusplus
}
#endif

#endif
