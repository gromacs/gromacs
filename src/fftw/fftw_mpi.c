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

#include <stdio.h>
#include <math.h>

#include <fftw_mpi.h>
#include <fftw-int.h>

/************************** Twiddle Factors *****************************/

/* To conserve space, we share twiddle factor arrays between forward and
   backward plans and plans of the same size (just as in the uniprocessor
   transforms). */

static fftw_mpi_twiddle *fftw_mpi_twiddles = NULL;

static fftw_mpi_twiddle *fftw_mpi_create_twiddle(int rows, int rowstart,
						 int cols, int n)
{
     fftw_mpi_twiddle *tw = fftw_mpi_twiddles;

     while (tw && (tw->rows != rows || tw->rowstart != rowstart ||
		   tw->cols != cols || tw->n != n))
	  tw = tw->next;

     if (tw) {
	  tw->refcount++;
	  return tw;
     }

     tw = (fftw_mpi_twiddle *) fftw_malloc(sizeof(fftw_mpi_twiddle));
     tw->rows = rows;
     tw->rowstart = rowstart;
     tw->cols = cols;
     tw->n = n;
     tw->refcount = 1;
     tw->next = fftw_mpi_twiddles;

     {
	  fftw_complex *W = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *
							 rows * (cols - 1));
	  int j, i;
	  FFTW_TRIG_REAL twoPiOverN = FFTW_K2PI / (FFTW_TRIG_REAL) n;

	  for (j = 0; j < rows; ++j)
	       for (i = 1; i < cols; ++i) {
		    int k = (j * (cols - 1) - 1) + i;
		    FFTW_TRIG_REAL
			 ij = (FFTW_TRIG_REAL) (i * (j + rowstart));
		    c_re(W[k]) = FFTW_TRIG_COS(twoPiOverN * ij);
		    c_im(W[k]) = FFTW_FORWARD * FFTW_TRIG_SIN(twoPiOverN * ij);
	       }
	  
	  tw->W = W;
     }

     fftw_mpi_twiddles = tw;

     return tw;
}

static void fftw_mpi_destroy_twiddle(fftw_mpi_twiddle *tw)
{
     if (tw) {
	  tw->refcount--;
	  if (tw->refcount == 0) {
	       /* delete tw from fftw_mpi_twiddles list: */
	       if (fftw_mpi_twiddles == tw)
		    fftw_mpi_twiddles = tw->next;
	       else {
		    fftw_mpi_twiddle *prev = fftw_mpi_twiddles;
		    
		    if (!prev)
			 fftw_mpi_die("unexpected empty MPI twiddle list");
		    while (prev->next && prev->next != tw)
			 prev = prev->next;
		    if (prev->next != tw)
			 fftw_mpi_die("tried to destroy unknown MPI twiddle");
		    prev->next = tw->next;
	       }
	       
	       fftw_free(tw->W);
	       fftw_free(tw);
	  }
     }
}

/* multiply the array in d (of size tw->cols * n_fields) by the row cur_row
   of the twiddle factors pointed to by tw, given the transform direction. */
static void fftw_mpi_mult_twiddles(fftw_complex *d, int n_fields,
				   int cur_row,
				   fftw_mpi_twiddle *tw,
				   fftw_direction dir)
{
     int cols = tw->cols;
     fftw_complex *W = tw->W + cur_row * (cols - 1);
     int j;

     if (dir == FFTW_FORWARD) {
          if (n_fields > 1)
               for (j = 1; j < cols; ++j) {
                    fftw_real
                         w_re = c_re(W[j-1]),
                         w_im = c_im(W[j-1]);
                    int f;

                    for (f = 0; f < n_fields; ++f) {
                         fftw_real
                              d_re = c_re(d[j*n_fields + f]),
                              d_im = c_im(d[j*n_fields + f]);
                         c_re(d[j*n_fields + f]) = w_re * d_re - w_im * d_im;
                         c_im(d[j*n_fields + f]) = w_re * d_im + w_im * d_re;
                    }
               }
          else
               for (j = 1; j < cols; ++j) {
                    fftw_real w_re = c_re(W[j-1]),
                         w_im = c_im(W[j-1]),
                         d_re = c_re(d[j]),
                         d_im = c_im(d[j]);
                    c_re(d[j]) = w_re * d_re - w_im * d_im;
                    c_im(d[j]) = w_re * d_im + w_im * d_re;
               }
     }
     else {  /* FFTW_BACKWARDS */
	  /* same as above, except that W is complex-conjugated: */
          if (n_fields > 1)
               for (j = 1; j < cols; ++j) {
                    fftw_real
                         w_re = c_re(W[j-1]),
                         w_im = c_im(W[j-1]);
                    int f;

                    for (f = 0; f < n_fields; ++f) {
                         fftw_real
                              d_re = c_re(d[j*n_fields + f]),
                              d_im = c_im(d[j*n_fields + f]);
                         c_re(d[j*n_fields + f]) = w_re * d_re + w_im * d_im;
                         c_im(d[j*n_fields + f]) = w_re * d_im - w_im * d_re;
                    }
               }
          else
               for (j = 1; j < cols; ++j) {
                    fftw_real w_re = c_re(W[j-1]),
                         w_im = c_im(W[j-1]),
                         d_re = c_re(d[j]),
                         d_im = c_im(d[j]);
                    c_re(d[j]) = w_re * d_re + w_im * d_im;
                    c_im(d[j]) = w_re * d_im - w_im * d_re;
               }
     }
}

/***************************** Plan Creation ****************************/

/* return the factor of n closest to sqrt(n): */
static int find_sqrt_factor(int n)
{
     int i = sqrt(n) + 0.5;
     int i2 = i - 1;
     
     while (i2 > 0) {
	  if (n % i2 == 0)
	       return i2;
	  if (n % i == 0)
	       return i;
	  ++i; --i2;
     }
     return 1; /* n <= 1 */
}

/* find the "best" r to divide n by for the FFT decomposition.  Ideally,
   we would like both r and n/r to be divisible by the number of 
   processes (for optimum load-balancing).  Also, pick r to be close
   to sqrt(n) if possible. */
static int find_best_r(int n, MPI_Comm comm)
{
     int n_pes;

     MPI_Comm_size(comm, &n_pes);

     if (n % n_pes == 0) {
	  n /= n_pes;
	  if (n % n_pes == 0)
	       return (n_pes * find_sqrt_factor(n / n_pes));
	  else
	       return (n_pes * find_sqrt_factor(n));
     }
     else
	  return find_sqrt_factor(n);
}

#define MAX2(a,b) ((a) > (b) ? (a) : (b))

fftw_mpi_plan fftw_mpi_create_plan(MPI_Comm comm,
				   int n, fftw_direction dir, int flags)
{
     fftw_mpi_plan p;
     int i, r, m;

     p = (fftw_mpi_plan) fftw_malloc(sizeof(struct fftw_mpi_plan_struct));

     i = find_best_r(n, comm);
     if (dir == FFTW_FORWARD)
	  m = n / (r = i);
     else
	  r = n / (m = i);

     p->n = n;
     p->r = r;
     p->m = m;

     flags |= FFTW_IN_PLACE;
     p->flags = flags;
     p->dir = dir;

     p->pr = fftw_create_plan(r, dir, flags);
     p->pm = fftw_create_plan(m, dir, flags);
	  
     p->p_transpose = transpose_mpi_create_plan(m, r, comm);
     p->p_transpose_inv = transpose_mpi_create_plan(r, m, comm);

     transpose_mpi_get_local_size(r,
				  p->p_transpose_inv->my_pe,
				  p->p_transpose_inv->n_pes,
				  &p->local_r,
				  &p->local_r_start);
     transpose_mpi_get_local_size(m,
				  p->p_transpose->my_pe,
				  p->p_transpose->n_pes,
				  &p->local_m,
				  &p->local_m_start);

     if (dir == FFTW_FORWARD)
	  p->tw = fftw_mpi_create_twiddle(p->local_r, p->local_r_start, m, n);
     else
	  p->tw = fftw_mpi_create_twiddle(p->local_m, p->local_m_start, r, n);

     p->fft_work = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) *
						MAX2(m, r));

     return p;
}

/********************* Getting Local Size ***********************/

void fftw_mpi_local_sizes(fftw_mpi_plan p,
			  int *local_n,
			  int *local_start,
			  int *local_n_after_transform,
			  int *local_start_after_transform,
			  int *total_local_size)
{
     if (p) {
	  if (p->flags & FFTW_SCRAMBLED_INPUT) {
	       *local_n = p->local_r * p->m;
	       *local_start = p->local_r_start * p->m;
	  }
	  else {
	       *local_n = p->local_m * p->r;
	       *local_start = p->local_m_start * p->r;
	  }
	  
	  if (p->flags & FFTW_SCRAMBLED_OUTPUT) {
	       *local_n_after_transform = p->local_m * p->r;
	       *local_start_after_transform = p->local_m_start * p->r;
	  }
	  else {
	       *local_n_after_transform = p->local_r * p->m;
	       *local_start_after_transform = p->local_r_start * p->m;
	  }

	  *total_local_size =
	       transpose_mpi_get_local_storage_size(p->p_transpose->nx,
						    p->p_transpose->ny,
						    p->p_transpose->my_pe,
						    p->p_transpose->n_pes);
     }
}

static void fftw_mpi_fprint_plan(FILE *f, fftw_mpi_plan p)
{
     fprintf(f, "mpi plan:\n");
     fprintf(f, "m = %d plan:\n", p->m);
     fftw_fprint_plan(f, p->pm);
     fprintf(f, "r = %d plan:\n", p->r);
     fftw_fprint_plan(f, p->pr);
}

void fftw_mpi_print_plan(fftw_mpi_plan p)
{
     fftw_mpi_fprint_plan(stdout, p);
}

/********************** Plan Destruction ************************/

void fftw_mpi_destroy_plan(fftw_mpi_plan p)
{
     if (p) {
	  fftw_destroy_plan(p->pr);
	  fftw_destroy_plan(p->pm);
	  transpose_mpi_destroy_plan(p->p_transpose);
	  transpose_mpi_destroy_plan(p->p_transpose_inv);
	  fftw_mpi_destroy_twiddle(p->tw);
	  fftw_free(p->fft_work);
	  fftw_free(p);
     }
}

/******************** Computing the Transform *******************/

void fftw_mpi(fftw_mpi_plan p, int n_fields,
	      fftw_complex *local_data, fftw_complex *work)
{
     int i;
     int el_size = (sizeof(fftw_complex) / sizeof(TRANSPOSE_EL_TYPE))
                   * n_fields;
     fftw_complex *fft_work;
     fftw_direction dir;
     fftw_mpi_twiddle *tw;

     if (n_fields < 1)
	  return;

     if (!(p->flags & FFTW_SCRAMBLED_INPUT))
	  transpose_mpi(p->p_transpose, el_size,
			(TRANSPOSE_EL_TYPE *) local_data,
			(TRANSPOSE_EL_TYPE *) work);

     tw = p->tw;
     dir = p->dir;
     fft_work = work ? work : p->fft_work;

     /* For forward plans, we multiply by the twiddle factors here,
	before the second transpose.  For backward plans, we multiply
	by the twiddle factors after the second transpose.  We do
	this so that forward and backward transforms can share the
	same twiddle factor array (noting that m and r are swapped
	for the two directions so that the local sizes will be compatible). */

     {
	  int rows = p->local_r, cols = p->m;
	  fftw_plan p_fft = p->pm;

	  if (dir == FFTW_FORWARD) {
	       for (i = 0; i < rows; ++i) {
		    fftw_complex *d = local_data + i * (cols * n_fields);
		    
		    fftw(p_fft, n_fields, d, n_fields, 1, fft_work, 1, 0);
		    fftw_mpi_mult_twiddles(d, n_fields, i, tw, FFTW_FORWARD);
	       }
	  }
	  else {
	       if (n_fields > 1)
		    for (i = 0; i < rows; ++i)
			 fftw(p_fft, n_fields, local_data + i*(cols*n_fields),
			      n_fields, 1, fft_work, 1, 0);
	       else
		    fftw(p_fft, rows, local_data, 1, cols, fft_work, 1, 0);
	  }
     }

     transpose_mpi(p->p_transpose_inv, el_size,
		   (TRANSPOSE_EL_TYPE *) local_data,
		   (TRANSPOSE_EL_TYPE *) work);

     {
	  int rows = p->local_m, cols = p->r;
	  fftw_plan p_fft = p->pr;

	  if (dir == FFTW_BACKWARD) {
	       for (i = 0; i < rows; ++i) {
		    fftw_complex *d = local_data + i * (cols * n_fields);
		    
		    fftw_mpi_mult_twiddles(d, n_fields, i, tw, FFTW_BACKWARD);
		    fftw(p_fft, n_fields, d, n_fields, 1, fft_work, 1, 0);
	       }
	  }
	  else {
	       if (n_fields > 1)
		    for (i = 0; i < rows; ++i)
			 fftw(p_fft, n_fields, local_data + i*(cols*n_fields),
			      n_fields, 1, fft_work, 1, 0);
	       else
		    fftw(p_fft, rows, local_data, 1, cols, fft_work, 1, 0);
	  }
     }

     if (!(p->flags & FFTW_SCRAMBLED_OUTPUT))
	  transpose_mpi(p->p_transpose, el_size,
			(TRANSPOSE_EL_TYPE *) local_data,
			(TRANSPOSE_EL_TYPE *) work);

     /* Yes, we really had to do three transposes...sigh. */
}


