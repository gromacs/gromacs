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

#include <stdlib.h>

#include <mpi.h>

#include "transpose_mpi.h"
#include "fftwnd_mpi.h"

/***************************** Plan Creation ****************************/

fftwnd_mpi_plan fftwnd_mpi_create_plan(MPI_Comm comm,
				       int rank, const int *n,
				       fftw_direction dir,
				       int flags)
{
    fftwnd_mpi_plan p;

    if (rank < 2)
	return 0;

    p = (fftwnd_mpi_plan) fftw_malloc(sizeof(fftwnd_mpi_aux_data));
    p->p_fft = 0;
    p->p_transpose = 0;
    p->p_transpose_inv = 0;

    p->p_fft = fftwnd_create_plan(rank, n, dir, flags | FFTW_IN_PLACE);
    if (!p->p_fft)
	fftwnd_mpi_destroy_plan(p);

    p->p_transpose = transpose_mpi_create_plan(n[0], n[1], comm);
    if (!p->p_transpose)
	fftwnd_mpi_destroy_plan(p);

    p->p_transpose_inv = transpose_mpi_create_plan(n[1], n[0], comm);
    if (!p->p_transpose_inv)
	fftwnd_mpi_destroy_plan(p);

    return p;
}

fftwnd_mpi_plan fftw2d_mpi_create_plan(MPI_Comm comm,
				       int nx, int ny,
				       fftw_direction dir, int flags)
{
    int n[2];

    n[0] = nx;
    n[1] = ny;

    return fftwnd_mpi_create_plan(comm, 2, n, dir, flags);
}

fftwnd_mpi_plan fftw3d_mpi_create_plan(MPI_Comm comm,
				       int nx, int ny, int nz,
				       fftw_direction dir, int flags)
{
    int n[3];

    n[0] = nx;
    n[1] = ny;
    n[2] = nz;

    return fftwnd_mpi_create_plan(comm, 3, n, dir, flags);
}

/********************** Plan Destruction ************************/

void fftwnd_mpi_destroy_plan(fftwnd_mpi_plan p)
{
    if (p) {
	if (p->p_fft)
	    fftwnd_destroy_plan(p->p_fft);
	if (p->p_transpose)
	    transpose_mpi_destroy_plan(p->p_transpose);
	if (p->p_transpose_inv)
	    transpose_mpi_destroy_plan(p->p_transpose_inv);
	fftw_free(p);
    }
}

/********************* Getting Local Size ***********************/

void fftwnd_mpi_local_sizes(fftwnd_mpi_plan p,
			    int *local_nx,
			    int *local_x_start,
			    int *local_ny_after_transpose,
			    int *local_y_start_after_transpose,
			    int *total_local_size)
{
    if (p) {
	transpose_mpi_get_local_size(p->p_transpose->nx,
				     p->p_transpose->my_pe,
				     p->p_transpose->n_pes,
				     local_nx,
				     local_x_start);
	transpose_mpi_get_local_size(p->p_transpose->ny,
				     p->p_transpose->my_pe,
				     p->p_transpose->n_pes,
				     local_ny_after_transpose,
				     local_y_start_after_transpose);
	*total_local_size =
	    transpose_mpi_get_local_storage_size(p->p_transpose->nx,
						 p->p_transpose->ny,
						 p->p_transpose->my_pe,
						 p->p_transpose->n_pes);

	*total_local_size *= p->p_fft->n_after[1];
    }
}

/******************** Computing the Transform *******************/

void fftwnd_mpi(fftwnd_mpi_plan p, int n_fields, fftw_complex * local_data,
		fftwnd_mpi_output_order output_order)
{
    int fft_iter;
    int j, i;
    fftw_plan *plans;
    int *n, *n_before, *n_after;
    int nb, nx, ny, local_nx, local_ny;
    fftw_complex *work;
    int rank;

    if (n_fields <= 0)
	return;

    plans = p->p_fft->plans;
    n = p->p_fft->n;
    n_before = p->p_fft->n_before;
    n_after = p->p_fft->n_after;
    work = p->p_fft->work;
    rank = p->p_fft->rank;
    nx = n[0];
    ny = n[1];
    local_nx = p->p_transpose->local_nx;
    local_ny = p->p_transpose->local_ny;

    for (fft_iter = 0; fft_iter < n_fields; ++fft_iter) {
	/* do last dimension: */
	nb = (n_before[rank - 1] / nx) * local_nx;
	fftw(plans[rank - 1], nb,
	     local_data + fft_iter, n_fields, n[rank - 1] * n_fields,
	     work, 1, 0);

	/* do other dimensions, except for first: */
	for (j = 1; j < rank - 1; ++j) {
	    nb = (n_before[j] / nx) * local_nx;
	    for (i = 0; i < nb; ++i)
		fftw(plans[j], n_after[j],
		local_data + fft_iter + i * n_fields * n[j] * n_after[j],
		     n_fields * n_after[j], n_fields,
		     work, 1, 0);
	}
    }

    /* transpose the first two indices: */
    transpose_mpi(p->p_transpose, (fftw_real *) local_data,
		  2 * n_fields * n_after[1]);

    for (fft_iter = 0; fft_iter < n_fields; ++fft_iter) {
	/* do first dimension: */
	nb = local_ny;
	for (i = 0; i < nb; ++i)
	    fftw(plans[0], n_after[1],
		 local_data + fft_iter + i * n_fields * n[0] * n_after[1],
		 n_fields * n_after[1], n_fields,
		 work, 1, 0);
    }

    /* transpose back, if desired: */
    if (output_order == FFTW_NORMAL_ORDER)
	transpose_mpi(p->p_transpose_inv, (fftw_real *) local_data,
		      2 * n_fields * n_after[1]);
}
