/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2005
 * David van der Spoel, Erik Lindahl, University of Groningen.
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

#ifndef _gmx_parallel_3dfft_h_
#define _gmx_parallel_3dfft_h_


#include "types/simple.h"
#include "types/commrec.h"
#include "gmxcomplex.h"
#include "gmx_fft.h"

typedef struct gmx_parallel_3dfft *
    gmx_parallel_3dfft_t;



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
 *  \param ngridx         Global number of grid cells in the x direction. Must be
 *                        divisible by the number of nodes.
 *  \param ngridy         Global number of grid cells in the y direction. Must be
 *                        divisible by the number of nodes.
 *  \param ngridz         Global number of grid cells in the z direction.
 *  \param node2slab      Node id to slab index array, can be NULL.
 *  \param slab2grid_x    Slab index to grid_x array (nnodes+1), can be NULL.
 *  \param comm           MPI communicator, must have been initialized.
 *  \param bReproducible  Try to avoid FFT timing optimizations and other stuff
 *                        that could make results differ for two runs with
 *                        identical input (reproducibility for debugging).
 *  \param nthreads       Run in parallel using n threads
 *
 *  \return 0 or a standard error code.
 */
int
    gmx_parallel_3dfft_init   (gmx_parallel_3dfft_t *    pfft_setup,
                               ivec                      ndata,
                               real **real_data,
                               t_complex **complex_data,
                               MPI_Comm                  comm[2],
                               int *                     slab2index_major,
                               int *                     slab2index_minor,
                               gmx_bool                  bReproducible,
                               int                       nthreads);





/*! \brief Get direct space grid index limits
 */
int
gmx_parallel_3dfft_real_limits(gmx_parallel_3dfft_t      pfft_setup,
                               ivec                      local_ndata,
                               ivec                      local_offset,
                               ivec                      local_size);


/*! \brief Get reciprocal space grid index limits
 */
int
gmx_parallel_3dfft_complex_limits(gmx_parallel_3dfft_t      pfft_setup,
                                  ivec                      complex_order,
                                  ivec                      local_ndata,
                                  ivec                      local_offset,
                                  ivec                      local_size);


int
gmx_parallel_3dfft_execute(gmx_parallel_3dfft_t    pfft_setup,
                           enum gmx_fft_direction  dir,
                           void *                  in_data,
                           void *                  out_data,
                           int                     thread,
                           gmx_wallcycle_t         wcycle);


/*! \brief Release all data in parallel fft setup
 *
 *  All temporary storage and FFT plans are released. The structure itself
 *  is not released, but the contents is invalid after this call.
 *
 *  \param pfft_setup Parallel 3dfft setup.
 *  \param in_data    Input data.
 *  \param out_data   Output data.
 *  \param thread     Thread index of the calling thread, i.e. index to the part
 *                    of the data operated on last by the calling thread. This
 *                    is needed to start the FFT without an OpenMP barrier.
 *  \param wcycle     Wall cycle counters.
 *
 *  \return 0 or a standard error code.
 */
int
gmx_parallel_3dfft_destroy(gmx_parallel_3dfft_t    pfft_setup);




#endif /* _gmx_parallel_3dfft_h_ */
