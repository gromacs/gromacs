/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 *  \brief Declare interface for GPU execution for NBNXN module
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_NBNXN_GPU_H
#define GMX_MDLIB_NBNXN_GPU_H

#include "gromacs/gmxlib/gpu_utils/gpu_macros.h"
#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/mdlib/nbnxn_gpu_types.h"

#ifdef __cplusplus
extern "C" {
#endif

struct nbnxn_atomdata_t;

/*! \brief
 * Launch asynchronously the nonbonded force calculations.
 *
 *  This consists of the following (async) steps launched:
 *  - upload x and q;
 *  - upload shift vector;
 *  - launch kernel;
 *  The local and non-local interaction calculations are launched in two
 *  separate streams.
 */
GPU_FUNC_QUALIFIER
void nbnxn_gpu_launch_kernel(gmx_nbnxn_gpu_t gmx_unused               *nb,
                             const struct nbnxn_atomdata_t gmx_unused *nbdata,
                             int gmx_unused                            flags,
                             int gmx_unused                            iloc) GPU_FUNC_TERM

/*! \brief
 * Launch asynchronously the download of nonbonded forces from the GPU
 * (and energies/shift forces if required).
 */
GPU_FUNC_QUALIFIER
void nbnxn_gpu_launch_cpyback(gmx_nbnxn_gpu_t  gmx_unused              *nb,
                              const struct nbnxn_atomdata_t gmx_unused *nbatom,
                              int                    gmx_unused         flags,
                              int                    gmx_unused         aloc) GPU_FUNC_TERM

/*! \brief
 * Wait for the asynchronously launched nonbonded calculations and data
 * transfers to finish.
 */
GPU_FUNC_QUALIFIER
void nbnxn_gpu_wait_for_gpu(gmx_nbnxn_gpu_t gmx_unused               *nb,
                            const struct nbnxn_atomdata_t gmx_unused *nbatom,
                            int                    gmx_unused         flags,
                            int                    gmx_unused         aloc,
                            real                   gmx_unused        *e_lj,
                            real                   gmx_unused        *e_el,
                            rvec                   gmx_unused        *fshift) GPU_FUNC_TERM

/*! \brief Selects the Ewald kernel type, analytical or tabulated, single or twin cut-off. */
GPU_FUNC_QUALIFIER
int nbnxn_gpu_pick_ewald_kernel_type(bool gmx_unused bTwinCut) GPU_FUNC_TERM_WITH_RETURN(-1)

#ifdef __cplusplus
}
#endif

#endif
