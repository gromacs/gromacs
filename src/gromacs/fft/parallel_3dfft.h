/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#ifndef GMX_FFT_PARALLEL_3DFFT_H
#define GMX_FFT_PARALLEL_3DFFT_H

#include "gromacs/fft/fft.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/real.h"

typedef struct gmx_parallel_3dfft* gmx_parallel_3dfft_t;


/*! \brief Initialize parallel MPI-based 3D-FFT.
 *
 *  This routine performs real-to-complex and complex-to-real parallel 3D FFTs,
 *  but not complex-to-complex.
 *
 *  The routine is optimized for small-to-medium size FFTs used for PME and
 *  PPPM algorithms, and do allocate extra workspace whenever it might improve
 *  performance.
 *
 *  \param pfft_setup     Pointer to parallel 3dfft setup structure, previously
 *                        allocated or with automatic storage.
 *  \param ndata          Number of grid cells in each direction
 *  \param real_data      Real data. Input for forward and output for backward.
 *  \param complex_data   Complex data.
 *  \param comm           MPI communicator for both parallelization axis.
 *                        Needs to be either initialized or MPI_NULL for
 *                        no parallelization in that axis.
 *  \param bReproducible  Try to avoid FFT timing optimizations and other stuff
 *                        that could make results differ for two runs with
 *                        identical input (reproducibility for debugging).
 *  \param nthreads       Run in parallel using n threads
 *  \param realGridAllocation  Whether to make real grid use allocation pinned for GPU transfers.
 *                             Only used in PME mixed CPU+GPU mode.
 *
 *  \return 0 or a standard error code.
 */
int gmx_parallel_3dfft_init(gmx_parallel_3dfft_t* pfft_setup,
                            const ivec            ndata,
                            real**                real_data,
                            t_complex**           complex_data,
                            MPI_Comm              comm[2],
                            gmx_bool              bReproducible,
                            int                   nthreads,
                            gmx::PinningPolicy realGridAllocation = gmx::PinningPolicy::CannotBePinned);


/*! \brief Get direct space grid index limits
 */
int gmx_parallel_3dfft_real_limits(gmx_parallel_3dfft_t pfft_setup,
                                   ivec                 local_ndata,
                                   ivec                 local_offset,
                                   ivec                 local_size);


/*! \brief Get reciprocal space grid index limits
 */
int gmx_parallel_3dfft_complex_limits(gmx_parallel_3dfft_t pfft_setup,
                                      ivec                 complex_order,
                                      ivec                 local_ndata,
                                      ivec                 local_offset,
                                      ivec                 local_size);


int gmx_parallel_3dfft_execute(gmx_parallel_3dfft_t   pfft_setup,
                               enum gmx_fft_direction dir,
                               int                    thread,
                               gmx_wallcycle*         wcycle);


/*! \brief Release all data in parallel fft setup
 *
 *  All temporary storage and FFT plans are released. The structure itself
 *  is not released, but the contents is invalid after this call.
 *
 *  \param pfft_setup Parallel 3dfft setup.
 *
 *  \return 0 or a standard error code.
 */
int gmx_parallel_3dfft_destroy(gmx_parallel_3dfft_t pfft_setup);

#endif
