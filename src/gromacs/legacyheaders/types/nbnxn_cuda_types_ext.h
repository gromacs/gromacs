/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
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

#ifndef NBNXN_CUDA_TYPES_EXT_H
#define NBNXN_CUDA_TYPES_EXT_H

#ifdef __cplusplus
extern "C" {
#endif

/* Abstract types */
/* CUDA nonbonded structure */
typedef struct nbnxn_cuda *nbnxn_cuda_ptr_t;
/* CUDA GPU device info */
typedef struct cuda_dev_info *cuda_dev_info_ptr_t;

/* Types defined for the structs below. */
typedef struct wallclock_gpu wallclock_gpu_t;
typedef struct nbnxn_cuda_ktime nbnxn_cuda_ktime_t;

/* Nonbonded kernel time and call count. */
struct nbnxn_cuda_ktime
{
    double  t;
    int     c;
};

/* GPU timings for kernels and H2d/D2H transfers. */
struct wallclock_gpu
{
    nbnxn_cuda_ktime_t ktime[2][2]; /* table containing the timings of the four
                                       version of the nonbonded kernels: force-only,
                                       force+energy, force+pruning, and force+energy+pruning */
    double  nb_h2d_t;               /* host to device transfer time in nb calculation  */
    double  nb_d2h_t;               /* device to host transfer time in nb calculation */
    int     nb_c;                   /* total call count of the nonbonded gpu operations */
    double  pl_h2d_t;               /* pair search step host to device transfer time */
    int     pl_h2d_c;               /* pair search step  host to device transfer call count */
};

#ifdef __cplusplus
}
#endif

#endif /* NBNXN_CUDA_TYPES_EXT_H */
