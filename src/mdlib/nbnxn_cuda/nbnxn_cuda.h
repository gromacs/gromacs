/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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

#ifndef NBNXN_CUDA_H
#define NBNXN_CUDA_H

#include "types/nbnxn_cuda_types_ext.h"

#ifdef GMX_GPU
#define FUNC_TERM ;
#else
#define FUNC_TERM {}
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*! Launch asynchronously the nonbonded force calculations.
 *  This consists of the following (async) steps launched:
 *  - upload x and q;
 *  - upload shift vector;
 *  - launch kernel;
 *  The local and non-local interaction calculations are launched in two
 *  separate streams.
 */
void nbnxn_cuda_launch_kernel(nbnxn_cuda_ptr_t        cu_nb,
                              const nbnxn_atomdata_t *nbdata,
                              int                     flags,
                              int                     iloc) FUNC_TERM

/*! Launch asynchronously the download of nonbonded forces from the GPU
 *  (and energies/shift forces if required).
 */
void nbnxn_cuda_launch_cpyback(nbnxn_cuda_ptr_t        cu_nb,
                               const nbnxn_atomdata_t *nbatom,
                               int                     flags,
                               int                     aloc) FUNC_TERM

/*! Wait for the asynchronously launched nonbonded calculations and data
 *  transfers to finish.
 */
void nbnxn_cuda_wait_gpu(nbnxn_cuda_ptr_t cu_nb,
                         const nbnxn_atomdata_t * nbatom,
                         int flags, int aloc,
                         real *e_lj, real *e_el,
                         rvec *fshift) FUNC_TERM

#ifdef __cplusplus
}
#endif

#undef FUNC_TERM

#endif /* NBNXN_CUDA_H */
