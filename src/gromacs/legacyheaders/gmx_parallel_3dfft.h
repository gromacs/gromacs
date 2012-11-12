/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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


#include "types/simple.h"
#include "types/commrec.h"
#include "gmxcomplex.h"
#include "gmx_fft.h"

#ifdef __cplusplus
extern "C" {
#endif

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
 *  \param ndata          Number of grid cells in each direction
 *  \param real_data      Real data. Input for forward and output for backward.
 *  \param complex_data   Complex data.
 *  \param comm           MPI communicator for both parallelization axis.
 *                        Needs to be either initialized or MPI_NULL for
 *                        no parallelization in that axis.
 *  \param slab2index_major Not used
 *  \param slab2index_minor Not used
 *  \param bReproducible  Try to avoid FFT timing optimizations and other stuff
 *                        that could make results differ for two runs with
 *                        identical input (reproducibility for debugging).
 *  \param nthreads       Run in parallel using n threads
 *
 *  \return 0 or a standard error code.
 */
int
    gmx_parallel_3dfft_init   (gmx_parallel_3dfft_t *    pfft_setup,
                               ivec ndata,
                               real **real_data,
                               t_complex **complex_data,
                               MPI_Comm comm[2],
                               int *                     slab2index_major,
                               int *                     slab2index_minor,
                               gmx_bool bReproducible,
                               int nthreads);





/*! \brief Get direct space grid index limits
 */
int
gmx_parallel_3dfft_real_limits(gmx_parallel_3dfft_t pfft_setup,
                               ivec                 local_ndata,
                               ivec                 local_offset,
                               ivec                 local_size);


/*! \brief Get reciprocal space grid index limits
 */
int
gmx_parallel_3dfft_complex_limits(gmx_parallel_3dfft_t pfft_setup,
                                  ivec                 complex_order,
                                  ivec                 local_ndata,
                                  ivec                 local_offset,
                                  ivec                 local_size);


int
gmx_parallel_3dfft_execute(gmx_parallel_3dfft_t   pfft_setup,
                           enum gmx_fft_direction dir,
                           void                 * in_data,
                           void                 * out_data,
                           int                    thread,
                           gmx_wallcycle_t        wcycle);


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
gmx_parallel_3dfft_destroy(gmx_parallel_3dfft_t pfft_setup);

#ifdef __cplusplus
}
#endif


#endif /* _gmx_parallel_3dfft_h_ */

