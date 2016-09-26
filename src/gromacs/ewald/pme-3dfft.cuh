/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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

/*! \libinternal \file
 *  \brief Defines the CUDA 3D-FFT functions for PME.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#ifndef PME3DFFT_CUH
#define PME3DFFT_CUH

#include "gmxpre.h"

#include "pme-internal.h"

typedef struct gmx_parallel_3dfft_gpu *gmx_parallel_3dfft_gpu_t; // TODO: refactor

/*! \brief \internal
 * Initializes the CUDA FFT structure for performing real-to-complex and complex-to-real 3D FFT
 * on a PME grid of a given size.
 *
 * \param[in] pfft_setup            The CUDA FFT structure to be initialized.
 * \param[in] ndata                 The PME real-space grid size.
 * \param[in] pme                   The PME data structure.
 */
void pme_gpu_init_3dfft(gmx_parallel_3dfft_gpu_t *pfft_setup,
                        ivec                      ndata,
                        const gmx_pme_t          *pme);

/*! \brief \internal
 * Destroys the CUDA FFT structure.
 *
 * \param     pfft_setup            The CUDA FFT structure to be destroyed.
 */
void pme_gpu_destroy_3dfft(const gmx_parallel_3dfft_gpu_t &pfft_setup);

/*! \brief \internal
 *
 * Returns the grid dimensions of the local complex PME grid.
 *
 * \param[in]    pfft_setup            The CUDA FFT structure to be examined.
 * \param[out]   local_ndata           The numbers of complex elements in the local grid.
 * \param[out]   local_offset          The offsets of the local grid.
 * \param[out]   local_size            The numbers of complex elements (with padding) in the local grid.
 */
void pme_gpu_get_3dfft_complex_limits(const gmx_parallel_3dfft_gpu_t pfft_setup,
                                      ivec                           local_ndata,
                                      ivec                           local_offset,
                                      ivec                           local_size);

#endif // PME3DFFT_CUH
