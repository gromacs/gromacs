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

#include <math.h>
#include <stdio.h>
#include <time.h>

#include <mpi.h>
#include <fftwnd_mpi.h>

/* Maximum array sizes in test gets multiplied by this.  Change it
   to a bigger number (e.g. 10) if you want to run bigger tests
   and have enough memory. */
#define SIZE_MULT 3

#ifndef RAND_MAX
#define RAND_MAX 32767
#endif

#define NUM_ITER 20

#define USE_RANDOM  0		/* use superior random() function instead of rand() */

#if USE_RANDOM
long random(void);
void srandom(unsigned int seed);

#define rand random
#define srand srandom
#define RANDOM_MAX 2147483647
#else
#define RANDOM_MAX RAND_MAX
#endif

extern void test_fft(int rank, int istride, int ostride, long max_size,
		     long num_iters, short in_place, fftw_direction dir);


int main(int argc, char **argv)
{
     MPI_Init(&argc,&argv);

     if (argc > 1)
	  srand(atoi(argv[1]));
     else
	  srand(42);

    /* Test forward transforms: */
     
     test_fft(1, 1, 1, 64 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(1, 1, 1, 64 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
     test_fft(2, 1, 1, 32 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(2, 1, 1, 32 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
     test_fft(3, 1, 1, 10 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(3, 1, 1, 10 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
     test_fft(4, 1, 1, 4 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(4, 1, 1, 4 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
     test_fft(5, 1, 1, 2 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(5, 1, 1, 2 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
     
     test_fft(1, 2, 3, 64 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(1, 2, 3, 64 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
     test_fft(2, 2, 3, 32 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(2, 2, 3, 32 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
     test_fft(3, 2, 3, 10 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(3, 2, 3, 10 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
     test_fft(4, 2, 3, 4 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(4, 2, 3, 4 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
     test_fft(5, 2, 3, 2 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(5, 2, 3, 2 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
     
     test_fft(1, 3, 2, 64 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(1, 3, 2, 64 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
     test_fft(2, 3, 2, 32 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(2, 3, 2, 32 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
     test_fft(3, 3, 2, 10 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(3, 3, 2, 10 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
     test_fft(4, 3, 2, 4 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
     test_fft(4, 3, 2, 4 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
    test_fft(5, 3, 2, 2 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
    test_fft(5, 3, 2, 2 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);

    test_fft(1, 3, 3, 64 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
    test_fft(1, 3, 3, 64 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
    test_fft(2, 3, 3, 32 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
    test_fft(2, 3, 3, 32 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
    test_fft(3, 3, 3, 10 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
    test_fft(3, 3, 3, 10 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
    test_fft(4, 3, 3, 4 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
    test_fft(4, 3, 3, 4 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
    test_fft(5, 3, 3, 2 * SIZE_MULT, NUM_ITER, 0, FFTW_FORWARD);
    test_fft(5, 3, 3, 2 * SIZE_MULT, NUM_ITER, 1, FFTW_FORWARD);
    
    /* Test backward transforms */
    
    test_fft(1, 1, 1, 64 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(1, 1, 1, 64 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(2, 1, 1, 32 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(2, 1, 1, 32 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(3, 1, 1, 10 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(3, 1, 1, 10 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(4, 1, 1, 4 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(4, 1, 1, 4 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(5, 1, 1, 2 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(5, 1, 1, 2 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);

    test_fft(1, 2, 3, 64 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(1, 2, 3, 64 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(2, 2, 3, 32 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(2, 2, 3, 32 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(3, 2, 3, 10 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(3, 2, 3, 10 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(4, 2, 3, 4 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(4, 2, 3, 4 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(5, 2, 3, 2 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(5, 2, 3, 2 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);

    test_fft(1, 3, 2, 64 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(1, 3, 2, 64 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(2, 3, 2, 32 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(2, 3, 2, 32 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(3, 3, 2, 10 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(3, 3, 2, 10 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(4, 3, 2, 4 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(4, 3, 2, 4 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(5, 3, 2, 2 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(5, 3, 2, 2 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);

    test_fft(1, 3, 3, 64 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(1, 3, 3, 64 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(2, 3, 3, 32 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(2, 3, 3, 32 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(3, 3, 3, 10 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(3, 3, 3, 10 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(4, 3, 3, 4 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(4, 3, 3, 4 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);
    test_fft(5, 3, 3, 2 * SIZE_MULT, NUM_ITER, 0, FFTW_BACKWARD);
    test_fft(5, 3, 3, 2 * SIZE_MULT, NUM_ITER, 1, FFTW_BACKWARD);

    MPI_Finalize();

    return 0;
}

#define MAX_RANK 10

#define CLOSE(a,b) (fabs((a) - (b))*2.0/(fabs(a)+fabs(b)+0.1) < 1e-3)

void test_fft(int rank, int istride, int ostride,
	      long max_size,
	      long num_iters, short in_place,
	      fftw_direction dir)
{
    int dims[MAX_RANK];
    long dims_reverse[MAX_RANK];
    fftw_real *a1, *a2, *a3, *a1b, *a2b;
    long i, j;
    int local_nx, local_x_start, local_ny, local_y_start, local_n;
    int howmany;
    long n_total, n_rest, iter;
    fftwnd_mpi_plan plan;
    fftwnd_plan seq_plan;
    int n_pes,my_pe;

    if (rank > MAX_RANK || rank < 2)
	return;

    if (!in_place)
	 return;
    else
	 ostride = istride;

    MPI_Comm_size(MPI_COMM_WORLD,&n_pes);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_pe);

    if (((max_size+n_pes-1)/n_pes)*(n_pes-1) >= max_size) {
	 if (my_pe == 0)
	      printf("*** max size too small for this many processes!\n");
	 return;
    }

    if (my_pe == 0)
    printf("TESTING RANK %d, stride %d/%d, sign %+d (in-place fftwnd)...\n",
	   rank, istride, ostride, (int) dir);

    for (iter = 0; iter < num_iters; ++iter) {
	 if (my_pe == 0) {
	      printf("   Testing ");
	      fflush(stdout);
	 }
	n_total = 1;
	for (i = 0; i < rank; ++i) {
	     do {
		  dims[i] = (rand() % max_size) + 1;
	     } while (i < 2 && ((dims[i]+n_pes-1)/n_pes)*(n_pes-1) >= dims[i]);
	     if (my_pe == 0) {
		  if (i)
		       printf("x%d", dims[i]);
		  else
		       printf("%d", dims[i]);
	     }
	    dims_reverse[rank - 1 - i] = dims[i];
	    n_total *= dims[i];
	}
	n_rest = n_total / dims[0];

	if (my_pe == 0) {
	     printf("...");
	     fflush(stdout);
	}

	seq_plan = fftwnd_create_plan(rank, dims, dir, FFTW_ESTIMATE |
				      FFTW_IN_PLACE);

	if ((rank != 2 && rank != 3) || (rand() & 1)) {
	     plan =
		  fftwnd_mpi_create_plan(MPI_COMM_WORLD,
					 rank, dims, dir, 
					 FFTW_ESTIMATE | FFTW_IN_PLACE);
	}
	/* for rank 2 and 3 arrays, use the optional interface 1/2 the
	 * time: */
	else if (rank == 2) {
	     if (my_pe == 0)
	    printf("(opt. interface)...");
		plan =
		    fftw2d_mpi_create_plan(MPI_COMM_WORLD,
					   dims[0], dims[1], dir, 
					   FFTW_ESTIMATE | FFTW_IN_PLACE);
	} else if (rank == 3) {
	     if (my_pe == 0)
	    printf("(opt. interface)...");
	    plan =
		 fftw3d_mpi_create_plan(MPI_COMM_WORLD,
					dims[0], dims[1], dims[2], 
					dir, FFTW_ESTIMATE | FFTW_IN_PLACE);
	}
	if (!plan) {
	    printf("\nError creating plan on process %d!\n",my_pe);
	    exit(1);
	}

	fftwnd_mpi_local_sizes(plan,&local_nx,&local_x_start,
			       &local_ny,&local_y_start,&local_n);

	a1 = (fftw_real*)fftw_malloc(n_total * 2 * sizeof(fftw_real));
	a2 = (fftw_real*)fftw_malloc(n_total * 2 * sizeof(fftw_real) * istride * 2);
	if (in_place) {
	    a3 = a2;
	    ostride = istride;
	} else
	    a3 = (fftw_real*)fftw_malloc(n_total * 2 * sizeof(fftw_real) * ostride);

	if (!a1 || !a2 || !a3) {
	    printf("\nERROR: Out of memory (need at least %ld bytes)\n", (n_total * 2 + n_total * 4) * sizeof(fftw_real));
	    exit(1);
	}
	for (i = 0; i < n_total; ++i) {
	    a1[2 * i] = (rand() - (RANDOM_MAX >> 1)) * 2.0 / RANDOM_MAX;
	    a1[2 * i + 1] = (rand() - (RANDOM_MAX >> 1)) * 2.0 /
		RANDOM_MAX;
	    for (j = 0; j < istride; ++j) {
		a2[2 * (i * istride + j)] = a1[2 * i];
		a2[2 * (i * istride + j) + 1] = a1[2 * i + 1];
	    }
	}

	howmany = istride;
	if (howmany > ostride)
	    howmany = ostride;

	if (my_pe == 0) {
	     printf("fftwnd_mpi...");
	     fflush(stdout);
	}

	/* FFT using MPI */

	fftwnd_mpi(plan,howmany,
		   (fftw_complex *) (a2 + 
				     2*istride*local_x_start*n_rest),
		   FFTW_NORMAL_ORDER);

	if (my_pe == 0) {
	     printf("fftwnd...");
	     fflush(stdout);
	}

	if (seq_plan->is_in_place != FFTW_IN_PLACE) {
	    printf("\n\nCorrupted plan!  (Invalid in-place flag.)\n");
	    exit(1);
	}
	/* FFT using sequential fftw */

	fftwnd(seq_plan, 1, (fftw_complex *) a1, 1, 1, 0, 0, 0);

	/* Check results: */

	a2b = a3 + 2*istride*local_x_start*n_rest;
	a1b = a1 + 2*local_x_start*n_rest;

	for (i = 0; i < local_nx*n_rest; ++i)
	    for (j = 0; j < howmany; ++j) {
		if (!CLOSE(a1b[2 * i], a2b[2 * (i * ostride + j)])) {
		    printf("error!\nFFT's different in real part at i=%ld: %g vs. %g\n",
			   i, a1b[2 * i], a2b[2 * (i * ostride + j)]);
		    exit(1);
		}
		if (!CLOSE(a1b[2 * i + 1], a2b[2 * (i * ostride +
						  j) + 1])) {
		    printf("error!\nFFT's different in imag. part at i=%ld: %g vs. %g\n",
			i, a1b[2 * i + 1], a2b[2 * (i * ostride + j) + 1]);
		    exit(1);
		}
	    }

	fftw_free(a1);
	fftw_free(a2);
	if (a3 != a2)
	    fftw_free(a3);

	fftwnd_mpi_destroy_plan(plan);
	fftwnd_destroy_plan(seq_plan);

	if (my_pe == 0) {
	     printf("okay!\n");
	     fflush(stdout);
	}

	fftw_check_memory_leaks();
    }
}
