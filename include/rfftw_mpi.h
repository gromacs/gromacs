/*
 * Copyright (c) 1997-1999 Massachusetts Institute of Technology
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

#ifndef RFFTW_MPI_H
#define RFFTW_MPI_H

#include <fftw_mpi.h>
#include <rfftw.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/***********************************************************************/

typedef struct {
     fftw_plan p_fft_x;  /* plan for first dimension */
     rfftwnd_plan p_fft;  /* plan for subsequent dimensions */
     transpose_mpi_plan p_transpose, p_transpose_inv;
     fftw_complex *work; /* extra workspace, if needed */
} rfftwnd_mpi_plan_data;

typedef rfftwnd_mpi_plan_data *rfftwnd_mpi_plan;

extern rfftwnd_mpi_plan rfftwnd_mpi_create_plan(MPI_Comm comm,
					      int rank, const int *n,
					      fftw_direction dir,
					      int flags);
extern rfftwnd_mpi_plan rfftw2d_mpi_create_plan(MPI_Comm comm,
					      int nx, int ny,
					  fftw_direction dir, int flags);
extern rfftwnd_mpi_plan rfftw3d_mpi_create_plan(MPI_Comm comm,
					      int nx, int ny, int nz,
					  fftw_direction dir, int flags);

extern void rfftwnd_mpi_destroy_plan(rfftwnd_mpi_plan p);

extern void rfftwnd_mpi_local_sizes(rfftwnd_mpi_plan p,
				   int *local_nx,
				   int *local_x_start,
				   int *local_ny_after_transpose,
				   int *local_y_start_after_transpose,
				   int *total_local_size);

extern void rfftwnd_mpi(rfftwnd_mpi_plan p,
		       int n_fields,
		       fftw_real *local_data, fftw_real *work,
		       fftwnd_mpi_output_order output_order);

/***********************************************************************/

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* RFFTW_MPI_H */
