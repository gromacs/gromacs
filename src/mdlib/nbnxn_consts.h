/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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


/* The number of pair-search sub-cells per super-cell, used for GPU */
#define GPU_NSUBCELL_Z 2
#define GPU_NSUBCELL_Y 2
#define GPU_NSUBCELL_X 2
#define GPU_NSUBCELL   (GPU_NSUBCELL_Z*GPU_NSUBCELL_Y*GPU_NSUBCELL_X)
/* In the non-bonded GPU kernel we operate on cluster-pairs, not cells.
 * The number of cluster in a super-cluster matches the number of sub-cells
 * in a pair-search cell, so we introduce a new name for the same value.
 */
#define NBNXN_GPU_NCLUSTER_PER_SUPERCLUSTER  GPU_NSUBCELL

/* With CPU kernels the i-cluster size is always 4 atoms.
 * With x86 SIMD the j-cluster size can be 2, 4 or 8, otherwise 4.
 */
#define NBNXN_CPU_CLUSTER_I_SIZE       4

#define NBNXN_CPU_CLUSTER_I_SIZE_2LOG  2

/* With GPU kernels the cluster size is 8 atoms */
#define NBNXN_GPU_CLUSTER_SIZE         8

/* With GPU kernels we group cluster pairs in 4 to optimize memory usage */
#define NBNXN_GPU_JGROUP_SIZE  4

/* To avoid NaN when excluded atoms are at zero distance, we add a small
 * number to r^2. NBNXN_AVOID_SING_R2_INC^-3 should fit in real.
 */
#ifndef GMX_DOUBLE
#define NBNXN_AVOID_SING_R2_INC  1.0e-12f
#else
/* The double prec. x86 SIMD kernels use a single prec. invsqrt, so > 1e-38 */
#define NBNXN_AVOID_SING_R2_INC  1.0e-36
#endif

/* Coulomb force table size chosen such that it fits along the non-bonded
   parameters in the texture cache. */
#define GPU_EWALD_COULOMB_FORCE_TABLE_SIZE 1536


/* Strides for x/f with xyz and xyzq coordinate (and charge) storage */
#define STRIDE_XYZ   3
#define STRIDE_XYZQ  4
/* Size of packs of x, y or z with SSE/AVX packed coords/forces */
#define PACK_X4      4
#define PACK_X8      8
/* Strides for a pack of 4 and 8 coordinates/forces */
#define STRIDE_P4    (DIM*PACK_X4)
#define STRIDE_P8    (DIM*PACK_X8)

/* Index of atom a into the SSE/AVX coordinate/force array */
#define X4_IND_A(a)  (STRIDE_P4*((a) >> 2) + ((a) & (PACK_X4 - 1)))
#define X8_IND_A(a)  (STRIDE_P8*((a) >> 3) + ((a) & (PACK_X8 - 1)))


#ifdef __cplusplus
}
#endif

#endif
