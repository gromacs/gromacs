/* This is not part of fftw-2.1.2, but a hack to be able to
 * use several threads on each node with mpi
 * Erik 990726
 */

#include <rfftw_mpi_threads.h>


/******************** Computing the Transform *******************/

static void first_dim_aux(int nthreads,rfftwnd_mpi_plan p,
			  int n_fields, fftw_real *local_data)
{
     int local_ny = p->p_transpose->local_ny;
     int nx = p->p_fft_x->n;
     fftw_complex *work_1d = p->work ? p->work : p->p_fft->work;
     
     n_fields *= p->p_fft->n_after[0]; /* dimensions after y 
					  no longer need be considered
					  separately from n_fields */
     /* single threaded transforms don't parallelize very well,
      * so instead we should try to divide transforms over threads
      */
     if (n_fields > 1) {
	  fftw_plan p_fft_x = p->p_fft_x;
	  int fft_iter;
	  for (fft_iter = 0; fft_iter < local_ny; ++fft_iter)
	      fftw_threads(nthreads,p_fft_x, n_fields,
			   ((fftw_complex *) local_data)
			   + (nx * n_fields) * fft_iter, n_fields, 1,
			   work_1d, 1, 0);
     }
     else 
	 fftw_threads(nthreads,p->p_fft_x, local_ny,
		      (fftw_complex *) local_data, 1, nx, work_1d, 1, 0);
}
			  
static void other_dims_aux(int nthreads,rfftwnd_mpi_plan p,
			  int n_fields, fftw_real *local_data)
{
     int local_nx = p->p_transpose->local_nx;
     int n_after_x = p->p_fft->n[0] * p->p_fft->n_after[0];
     
     if (n_fields > 1) {
	  rfftwnd_plan p_fft = p->p_fft;
	  int fft_iter;
	  if (p_fft->dir == FFTW_REAL_TO_COMPLEX)
	       for (fft_iter = 0; fft_iter < local_nx; ++fft_iter)
		    rfftwnd_threads_real_to_complex(nthreads,p_fft, n_fields,
				 local_data
				 + (2 * n_after_x * n_fields) * fft_iter,
				 n_fields, 1,
				 NULL, 0, 0);
	  else
	       for (fft_iter = 0; fft_iter < local_nx; ++fft_iter)
		    rfftwnd_threads_complex_to_real(nthreads,p_fft, n_fields,
				 ((fftw_complex *) local_data)
				 + (n_after_x * n_fields) * fft_iter,
				 n_fields, 1,
				 NULL, 0, 0);
     }
     else {
	  if (p->p_fft->dir == FFTW_REAL_TO_COMPLEX)
	       rfftwnd_threads_real_to_complex(nthreads,p->p_fft, local_nx,
				       local_data, 1, 2*n_after_x,
				       NULL, 0, 0);
	  else
	       rfftwnd_threads_complex_to_real(nthreads,p->p_fft, local_nx,
				       (fftw_complex *) local_data,
				       1, n_after_x,
				       NULL, 0, 0);
     }
}
			  

void rfftwnd_mpi_threads(int nthreads,rfftwnd_mpi_plan p,
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
	  other_dims_aux(nthreads,p, n_fields, local_data);
	  
	  /* Second, transpose the first dimension with the second dimension
	     to bring the x dimension local to this process: */
	  transpose_mpi(p->p_transpose, el_size, 
			(TRANSPOSE_EL_TYPE *) local_data,
			(TRANSPOSE_EL_TYPE *) work);
	  
	  /* Third, transform the x dimension, which is now 
	     local and contiguous: */
	  first_dim_aux(nthreads,p, n_fields, local_data);
	  
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
	  first_dim_aux(nthreads,p, n_fields, local_data);
	  
	  /* Third, transpose the first dimension with the second dimension
	     to bring the others dimensions local to this process: */
	  transpose_mpi(p->p_transpose_inv, el_size, 
			(TRANSPOSE_EL_TYPE *) local_data,
			(TRANSPOSE_EL_TYPE *) work);

	  /* last, transform dimensions after the first, which are
	     local to this process: */
	  other_dims_aux(nthreads,p, n_fields, local_data);
     }
}
