/*
 * Copyright (c) 1997,1998 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef FFTWND_MPI_H
#define FFTWND_MPI_H

#include "fftw.h"
#include "transpose_mpi.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef struct {
    fftwnd_plan p_fft;
    transpose_mpi_plan p_transpose, p_transpose_inv;
} fftwnd_mpi_aux_data;

typedef fftwnd_mpi_aux_data *fftwnd_mpi_plan;

typedef enum {
    FFTW_NORMAL_ORDER,
    FFTW_TRANSPOSED_ORDER
} fftwnd_mpi_output_order;

extern fftwnd_mpi_plan fftwnd_mpi_create_plan(MPI_Comm comm,
					      int rank, const int *n,
					      fftw_direction dir,
					      int flags);
extern fftwnd_mpi_plan fftw2d_mpi_create_plan(MPI_Comm comm,
					      int nx, int ny,
					      fftw_direction dir, int flags);
extern fftwnd_mpi_plan fftw3d_mpi_create_plan(MPI_Comm comm,
					      int nx, int ny, int nz,
					      fftw_direction dir, int flags);

extern void fftwnd_mpi_destroy_plan(fftwnd_mpi_plan p);

extern void fftwnd_mpi_local_sizes(fftwnd_mpi_plan p,
				   int *local_nx,
				   int *local_x_start,
				   int *local_ny_after_transpose,
				   int *local_y_start_after_transpose,
				   int *total_local_size);

extern void fftwnd_mpi(fftwnd_mpi_plan p,
		       int n_fields, fftw_complex * local_data,
		       fftwnd_mpi_output_order output_order);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* FFTWND_MPI_H */
