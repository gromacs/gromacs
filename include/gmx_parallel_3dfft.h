/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * $Id$
 *
 * Gromacs                               Copyright (c) 1991-2005
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#ifndef _gmx_parallel_3dfft_h_
#define _gmx_parallel_3dfft_h_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_MPI

#include "types/simple.h"
#include "gmxcomplex.h"
#include "gmx_fft.h"

/* We NEED MPI here. */
#include <mpi.h>

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
 *    
 *  \return 0 or a standard error code.
 */
int
gmx_parallel_3dfft_init   (gmx_parallel_3dfft_t *    pfft_setup,
                           int                       ngridx,
                           int                       ngridy,
                           int                       ngridz,
						   int                       *node2slab,
						   int                       *slab2grid_x,
                           MPI_Comm                  comm,
                           bool                      bReproducible);
                           




/*! \brief Get direct space grid index limits
 *
 *  The z dimension is never distributed. In the direct space, the x dimension
 *  is distributed over nodes, and after the real-to-complex FFT we work with
 *  a transposed grid where the y dimension is partitioned over nodes.
 *
 *  The node2slab array translates to node ids to slab indices,
 *  when NULL the slab ids are assumed to be identical to the node ids
 *  in the communicator comm.
 */
int
gmx_parallel_3dfft_limits(gmx_parallel_3dfft_t      pfft_setup,
                          int *                     local_x_start,
                          int *                     local_nx,
                          int *                     local_y_start,
                          int *                     local_ny);


int
gmx_parallel_transpose(t_complex *   data,
                       t_complex *   work,
                       int           nx,
                       int           ny,
                       int           local_x_start,
                       int           local_nx,
                       int           local_y_start,
                       int           local_ny,
                       int           nelem,
					   int           nnodes,
					   int           *node2slab,
                       MPI_Comm      comm);


/*! \brief Perform forward parallel MPI FFT.
 *
 *  Direction is either GMX_FFT_REAL_TO_COMPLEX or GMX_FFT_COMPLEX_TO_REAL.
 *
 *  If input and output arrays are separate there is no packing to consider.
 *  Input is simply nx*ny*nz in real, and output ny*nx*nzc in complex.
 *
 *  In they are identical we need to make sure there is room for the complex
 *  (length nzc=nz/2+1) in the array, so the _real_ space dimensions is
 *  always padded to nzc*2.
 *  In this case, the real dimensions are nx*ny*(nzc*2) while the complex
 *  dimensions is ny*nx*nzc (of type complex).
 *
 *  Note that the X and Y dimensions are transposed in the reciprocal space
 *  to avoid extra communication!
 *
 *  The node2slab array translates to node ids to slab indices,
 *  when NULL the slab ids are assumed to be identical to the node ids
 *  in the communicator comm.
 */
int
gmx_parallel_3dfft(gmx_parallel_3dfft_t    pfft_setup,
                   enum gmx_fft_direction  dir,
                   void *                  in_data,
                   void *                  out_data);



/*! \brief Release all data in parallel fft setup
 *
 *  All temporary storage and FFT plans are released. The structure itself
 *  is not released, but the contents is invalid after this call.
 *
 *  \param pfft_setup Parallel 3dfft setup.
 *
 *  \return 0 or a standard error code.
 */
int
gmx_parallel_3dfft_destroy(gmx_parallel_3dfft_t    pfft_setup);

#endif /* GMX_MPI */

#endif /* _gmx_parallel_3dfft_h_ */

