/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_NBNXN_CUDA_GPUBONDEDCUDA_H
#define GMX_MDLIB_NBNXN_CUDA_GPUBONDEDCUDA_H

#include "gromacs/mdlib/nbnxn_gpu_types.h"
#include "gromacs/ewald/pme-internal.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme-grid.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/enerdata.h"

#include "gromacs/mdlib/nbnxn_atomdata.h"

#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/pbcutil/mshift.h"

//CUDA threads per block
#define TPB_BONDED 256

//Macro for checking cuda errors following a cuda launch or api call
#define cudaCheckError() {                                                            \
    cudaError_t e=cudaGetLastError();                                                 \
    if(e!=cudaSuccess)                                                                \
    {                                                                                 \
        printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e)); \
        exit(0);                                                                      \
    }                                                                                 \
}
  
#if defined(__CUDACC__)
#include "nbnxn_cuda.h"
#include "nbnxn_cuda_types.h"
#include "cuda.h"
#include "cuda_runtime.h"
#endif

// new stuff for bonded
void update_gpu_bonded(const t_idef *idef,  const t_forcerec *fr, const matrix box,
                       const int size,  const t_mdatoms *md, gmx_grppairener_t *grppener); 
void do_bonded_gpu(t_forcerec *fr, const t_inputrec *ir, const t_idef *idef, 
                   int flags, const t_graph *graph , int natoms, rvec x[]);

void do_bonded_gpu_finalize(t_forcerec *fr, int flags, int natoms,
                            rvec *input_force, gmx_enerdata_t *enerd);

void reset_gpu_bonded(const int size, const int nener);

#endif



