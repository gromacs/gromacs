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

#include <rfftw_mpi.h>

/***************************** Plan Creation ****************************/

rfftwnd_mpi_plan rfftwnd_mpi_create_plan(MPI_Comm comm,
				       int rank, const int *n,
				       fftw_direction dir,
				       int flags)
{
    rfftwnd_mpi_plan p;

    if (rank < 2)
	return 0;

    p = (rfftwnd_mpi_plan) fftw_malloc(sizeof(rfftwnd_mpi_plan_data));
    p->p_fft_x = 0;
    p->p_fft = 0;
    p->p_transpose = 0;
    p->p_transpose_inv = 0;
    p->work = 0;

    p->p_fft_x = fftw_create_plan(n[0], dir, flags | FFTW_IN_PLACE);

    p->p_fft = rfftwnd_create_plan(rank-1, n+1, dir, flags | FFTW_IN_PLACE);
    if (!p->p_fft)
	rfftwnd_mpi_destroy_plan(p);

    p->p_transpose = transpose_mpi_create_plan(n[0], p->p_fft->n[0], comm);
    if (!p->p_transpose)
	rfftwnd_mpi_destroy_plan(p);

    p->p_transpose_inv = transpose_mpi_create_plan(p->p_fft->n[0], n[0], comm);
    if (!p->p_transpose_inv)
	rfftwnd_mpi_destroy_plan(p);

    if (n[0] > p->p_fft->nwork)
	 p->work = (fftw_complex *) fftw_malloc(n[0] * sizeof(fftw_complex));

    return p;
}

rfftwnd_mpi_plan rfftw2d_mpi_create_plan(MPI_Comm comm,
				         int nx, int ny,
				         fftw_direction dir, int flags)
{
    int n[2];

    n[0] = nx;
    n[1] = ny;

    return rfftwnd_mpi_create_plan(comm, 2, n, dir, flags);
}

rfftwnd_mpi_plan rfftw3d_mpi_create_plan(MPI_Comm comm,
			  	         int nx, int ny, int nz,
				         fftw_direction dir, int flags)
{
    int n[3];

    n[0] = nx;
    n[1] = ny;
    n[2] = nz;

    return rfftwnd_mpi_create_plan(comm, 3, n, dir, flags);
}

/********************** Plan Destruction ************************/

void rfftwnd_mpi_destroy_plan(rfftwnd_mpi_plan p)
{
    if (p) {
	if (p->p_fft_x)
	    fftw_destroy_plan(p->p_fft_x);
	if (p->p_fft)
	    rfftwnd_destroy_plan(p->p_fft);
	if (p->p_transpose)
	    transpose_mpi_destroy_plan(p->p_transpose);
	if (p->p_transpose_inv)
	    transpose_mpi_destroy_plan(p->p_transpose_inv);
	if (p->work)
	     fftw_free(p->work);
	fftw_free(p);
    }
}

/********************* Getting Local Size ***********************/

void rfftwnd_mpi_local_sizes(rfftwnd_mpi_plan p,
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

	*total_local_size *= p->p_fft->n_after[0];
	*total_local_size *= 2; /* return size in fftw_real's */

	if (p->p_fft->rank == 1 && p->p_fft->dir == FFTW_COMPLEX_TO_REAL) {
	     *local_ny_after_transpose *= 2;
	     *local_y_start_after_transpose *= 2;
	}
    }
}

/******************** Computing the Transform *******************/

static void first_dim_aux(rfftwnd_mpi_plan p,
			  int n_fields, fftw_real *local_data)
{
     int local_ny = p->p_transpose->local_ny;
     int nx = p->p_fft_x->n;
     fftw_complex *work_1d = p->work ? p->work : p->p_fft->work;
     
     n_fields *= p->p_fft->n_after[0]; /* dimensions after y 
					  no longer need be considered
					  separately from n_fields */
     if (n_fields > 1) {
	  fftw_plan p_fft_x = p->p_fft_x;
	  int fft_iter;
	  for (fft_iter = 0; fft_iter < local_ny; ++fft_iter)
	       fftw(p_fft_x, n_fields,
		    ((fftw_complex *) local_data)
		    + (nx * n_fields) * fft_iter, n_fields, 1,
		    work_1d, 1, 0);
     }
     else
	  fftw(p->p_fft_x, local_ny,
	       (fftw_complex *) local_data, 1, nx, work_1d, 1, 0);
}
			  
static void other_dims_aux(rfftwnd_mpi_plan p,
			  int n_fields, fftw_real *local_data)
{
     int local_nx = p->p_transpose->local_nx;
     int n_after_x = p->p_fft->n[0] * p->p_fft->n_after[0];
     
     if (n_fields > 1) {
	  rfftwnd_plan p_fft = p->p_fft;
	  int fft_iter;
	  if (p_fft->dir == FFTW_REAL_TO_COMPLEX)
	       for (fft_iter = 0; fft_iter < local_nx; ++fft_iter)
		    rfftwnd_real_to_complex(p_fft, n_fields,
				 local_data
				 + (2 * n_after_x * n_fields) * fft_iter,
				 n_fields, 1,
				 NULL, 0, 0);
	  else
	       for (fft_iter = 0; fft_iter < local_nx; ++fft_iter)
		    rfftwnd_complex_to_real(p_fft, n_fields,
				 ((fftw_complex *) local_data)
				 + (n_after_x * n_fields) * fft_iter,
				 n_fields, 1,
				 NULL, 0, 0);
     }
     else {
	  if (p->p_fft->dir == FFTW_REAL_TO_COMPLEX)
	       rfftwnd_real_to_complex(p->p_fft, local_nx,
				       local_data, 1, 2*n_after_x,
				       NULL, 0, 0);
	  else
	       rfftwnd_complex_to_real(p->p_fft, local_nx,
				       (fftw_complex *) local_data,
				       1, n_after_x,
				       NULL, 0, 0);
     }
}
			  

void rfftwnd_mpi(rfftwnd_mpi_plan p,
		 int n_fields, fftw_real *local_data, fftw_real *work,
		 fftwnd_mpi_output_order output_order)
{
     int el_size = (sizeof(fftw_complex) / sizeof(TRANSPOSE_EL_TYPE))
	           * n_fields * p->p_fft->n_after[0];
     
     if (n_fields <= 0)
	  return;

     if (p->p_fft->dir == FFTW_REAL_TO_COMPLEX) {
	  /* First, transform dimensions after the first, which are
	     local to this process: */	  
	  other_dims_aux(p, n_fields, local_data);
	  
	  /* Second, transpose the first dimension with the second dimension
	     to bring the x dimension local to this process: */
	  transpose_mpi(p->p_transpose, el_size, 
			(TRANSPOSE_EL_TYPE *) local_data,
			(TRANSPOSE_EL_TYPE *) work);
	  
	  /* Third, transform the x dimension, which is now 
	     local and contiguous: */
	  first_dim_aux(p, n_fields, local_data);
	  
	  /* transpose back, if desired: */
	  if (output_order == FFTW_NORMAL_ORDER)
	       transpose_mpi(p->p_transpose_inv, el_size,
			     (TRANSPOSE_EL_TYPE *) local_data,
			     (TRANSPOSE_EL_TYPE *) work);
     }
     else {  /* we have to do the steps in reverse order for c2r transform: */

	  /* NOTE: we assume that the same output_order is used for both
	     the forward and backward transforms: */

	  /* First, if necessary, transpose to get x dimension local: */
	  if (output_order == FFTW_NORMAL_ORDER)
	       transpose_mpi(p->p_transpose, el_size,
			     (TRANSPOSE_EL_TYPE *) local_data,
			     (TRANSPOSE_EL_TYPE *) work);
	  
	  /* Second, transform the x dimension, which is now 
	     local and contiguous: */
	  first_dim_aux(p, n_fields, local_data);
	  
	  /* Third, transpose the first dimension with the second dimension
	     to bring the others dimensions local to this process: */
	  transpose_mpi(p->p_transpose_inv, el_size, 
			(TRANSPOSE_EL_TYPE *) local_data,
			(TRANSPOSE_EL_TYPE *) work);

	  /* last, transform dimensions after the first, which are
	     local to this process: */
	  other_dims_aux(p, n_fields, local_data);
     }
}
