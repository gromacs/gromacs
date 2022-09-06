/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
#include "gromacs/nbnxm/simd_coulomb_functions.h"
#include "gromacs/nbnxm/simd_diagonal_masker.h"
#include "gromacs/nbnxm/simd_lennardjones_functions.h"
#include "gromacs/nbnxm/simd_load_store_functions.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/simd/vector_operations.h"
#ifdef CALC_COUL_EWALD
#    include "gromacs/math/utilities.h"
#endif

#include "config.h"

#include <cstdint>

#if !GMX_SIMD_HAVE_HSIMD_UTIL_REAL
#    error "Half-simd utility operations are required for the 2xNN kernels"
#endif

#ifndef GMX_SIMD_J_UNROLL_SIZE
#    error "Need to define GMX_SIMD_J_UNROLL_SIZE before including the 2xnn kernel common header file"
#endif

#define UNROLLI 4
#define UNROLLJ (GMX_SIMD_REAL_WIDTH / GMX_SIMD_J_UNROLL_SIZE)

static_assert(UNROLLI == c_nbnxnCpuIClusterSize, "UNROLLI should match the i-cluster size");

/* The stride of all the atom data arrays is equal to half the SIMD width */
#define STRIDE UNROLLJ

#if !defined GMX_NBNXN_SIMD_2XNN && !defined GMX_NBNXN_SIMD_4XN
#    error "Must define an NBNxN kernel flavour before including NBNxN kernel utility functions"
#endif

/* All functionality defines are set here, except for:
 * CALC_ENERGIES, ENERGY_GROUPS which are defined before.
 * CHECK_EXCLS, which is set just before including the inner loop contents.
 * The combination rule defines, LJ_COMB_GEOM or LJ_COMB_LB are currently
 * set before calling the kernel function. We might want to move that
 * to inside the n-loop and have a different combination rule for different
 * ci's, as no combination rule gives a 50% performance hit for LJ.
 */

/* We always calculate shift forces, because it's cheap anyhow */
#define CALC_SHIFTFORCES

/* Assumes all LJ parameters are identical */
/* #define FIX_LJ_C */
