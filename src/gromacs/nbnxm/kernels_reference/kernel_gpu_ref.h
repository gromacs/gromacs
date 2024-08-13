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
/*! \internal \file
 *
 * \brief
 * Declares GPU reference kernel
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_KERNELS_REFERENCE_KERNEL_GPU_REF_H
#define GMX_NBNXM_KERNELS_REFERENCE_KERNEL_GPU_REF_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

struct interaction_const_t;
struct t_forcerec;

namespace gmx
{
struct NbnxnPairlistGpu;
struct nbnxn_atomdata_t;
class StepWorkload;

//! Reference (slow) kernel for nb n vs n GPU type pair lists
void nbnxn_kernel_gpu_ref(const NbnxnPairlistGpu*    nbl,
                          const nbnxn_atomdata_t*    nbat,
                          const interaction_const_t* iconst,
                          ArrayRef<const RVec>       shiftvec,
                          const StepWorkload&        stepWork,
                          int                        clearF,
                          ArrayRef<real>             f,
                          real*                      fshift,
                          real*                      Vc,
                          real*                      Vvdw);

} // namespace gmx

#endif
