/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018, by the GROMACS development team, led by
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
 * \brief This file defines the PME CUDA-specific kernel parameter data structure.
 * \todo Rename the file (pme-gpu-types.cuh?), reconsider inheritance approach.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#ifndef GMX_EWALD_PME_CUH
#define GMX_EWALD_PME_CUH

#include "gromacs/math/vectypes.h" // for DIM

#include "pme-gpu-constants.h"
#include "pme-gpu-internal.h"      // for GridOrdering
#include "pme-gpu-types.h"

/*! \brief \internal
 * A single structure encompassing all the PME data used in CUDA kernels.
 * This inherits from PmeGpuKernelParamsBase and adds a couple cudaTextureObject_t handles,
 * which we would like to avoid in plain C++.
 */
struct PmeGpuCudaKernelParams : PmeGpuKernelParamsBase
{
    /* These are CUDA texture objects, related to the grid size. */
    /*! \brief CUDA texture object for accessing grid.d_fractShiftsTable */
    cudaTextureObject_t fractShiftsTableTexture;
    /*! \brief CUDA texture object for accessing grid.d_gridlineIndicesTable */
    cudaTextureObject_t gridlineIndicesTableTexture;
};

#endif
