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

#include <stdio.h>
#include <math.h>

#include <mpi.h>
#include <fftwnd_mpi.h>
#include <fftw-int.h>

#define NUM_ITER 90000000L

#define N_TESTS_3D 8

extern void initialize_fft_data(fftw_complex * arr, long n);

int main(int argc, char **argv)
{
    int n3[N_TESTS_3D][3] = 
    {
	 { 16, 16, 16 },
	 { 24, 24, 24 },
	 { 32, 32, 32 },
	 { 48, 48, 48 },
	 { 64, 64, 64 },
	 { 80, 80, 80 },
	 {100,100,100 },
	 { 128, 128, 128 },
/*         {256, 256, 256 }, */
    };
    int i, test, iter, max_iter;
    fftw_time start_t, end_t, init_t;
    fftw_complex *cin, *out;
    double time_scale,time1,time2;
    int max_size;
    fftwnd_plan plan_nd;
    fftwnd_mpi_plan plan_nd_mpi;
    int my_pe;

     MPI_Init(&argc,&argv);

     MPI_Comm_rank(MPI_COMM_WORLD,&my_pe);

    /*************** Benchmark fftwnd_mpi ****************/
     if (my_pe == 0)
    printf("\n");

    max_size = 0;
    for (i = 0; i < N_TESTS_3D; ++i)
	if (n3[i][0]*n3[i][1]*n3[i][2] > max_size)
	    max_size = n3[i][0]*n3[i][1]*n3[i][2];

    if (my_pe == 0) {
	 cin = (fftw_complex*)fftw_malloc(max_size * sizeof(fftw_complex));
	 
	 if (!cin) {
	      printf("Not enough memory!  At least %d bytes needed.\n",
		     max_size * sizeof(fftw_complex));
	      exit(1);
	 }
    }

     if (my_pe == 0)
    printf("%11s%12s%13s%13s%17s%13s\n", 
	   "Array Size", "FFTWND", "FFTWND_MPI",
	   "Speedup","FFTWND_MPI (T)","Speedup (T)");

    for (test = 0; test < N_TESTS_3D; ++test) {
	int N;
	int local_nx, local_x_start, local_ny, local_y_start, local_N;

	if (my_pe == 0)
	     plan_nd = fftwnd_create_plan(3,n3[test], FFTW_FORWARD, 
					  FFTW_IN_PLACE | FFTW_MEASURE);
	plan_nd_mpi = fftwnd_mpi_create_plan(MPI_COMM_WORLD,
					     3,n3[test], FFTW_FORWARD, 
					     FFTW_IN_PLACE | FFTW_MEASURE);

	fftwnd_mpi_local_sizes(plan_nd_mpi,&local_nx,&local_x_start,
                               &local_ny,&local_y_start,&local_N);

	if (my_pe != 0)
	     cin = (fftw_complex*)fftw_malloc(local_N * sizeof(fftw_complex));

	N = n3[test][0]*n3[test][1]*n3[test][2];
	local_N = local_nx * n3[test][1]*n3[test][2];

	max_iter = NUM_ITER / (N * log(2.0 * N));

	if (max_iter < 1)
	     max_iter = 1;

	time_scale = 1.0e6 / (max_iter * (log(N)/log(2.0) * N));

	if (my_pe == 0) {
	     initialize_fft_data(cin, N);
	     start_t = fftw_get_time();
	     for (iter = 0; iter < max_iter; ++iter)
		  initialize_fft_data(cin, N);
	     end_t = fftw_get_time();
	     init_t = fftw_time_diff(end_t,start_t);
	     
	     {
		  char s[20];
		  sprintf(s,"%dx%dx%d",n3[test][0],n3[test][1],n3[test][2]);
		  if (my_pe == 0)
		       printf("%11s",s);
	     }
	     fflush(stdout);
	     
	     /* Time FFTW: */
	     
	     initialize_fft_data(cin, N);
	     fftwnd(plan_nd, 1, cin, 1, 0, out, 1, 0);
	     start_t = fftw_get_time();
	     for (iter = 0; iter < max_iter; ++iter) {
		  initialize_fft_data(cin, N);
		  fftwnd(plan_nd, 1, cin, 1, 0, out, 1, 0);
	     }
	     end_t = fftw_get_time();
	     if (my_pe == 0) {
		  printf("%12g", 
			 time1=fftw_time_to_sec(fftw_time_diff(fftw_time_diff(end_t,start_t),init_t)) *
			 time_scale);
		  fflush(stdout);
	     }
	}
	else
	     time1 = 0;

	/* Time MPI FFTW: */

	initialize_fft_data(cin, local_N);
	start_t = fftw_get_time();
	for (iter = 0; iter < max_iter; ++iter)
	    initialize_fft_data(cin, local_N);
	end_t = fftw_get_time();
	init_t = fftw_time_diff(end_t,start_t);

	initialize_fft_data(cin, local_N);
	fftwnd_mpi(plan_nd_mpi, 1, cin, FFTW_NORMAL_ORDER);
	MPI_Barrier(MPI_COMM_WORLD);
	start_t = fftw_get_time();
	for (iter = 0; iter < max_iter; ++iter) {
	    initialize_fft_data(cin, local_N);
	    fftwnd_mpi(plan_nd_mpi, 1, cin, FFTW_NORMAL_ORDER);
	}
	end_t = fftw_get_time();
	if (my_pe == 0) {
	printf("%13g", time2=fftw_time_to_sec(fftw_time_diff(fftw_time_diff(end_t,start_t),init_t)) *
                       time_scale);
	printf("%13g",time1/time2);
	fflush(stdout);
	}

	/* Time MPI FFTW (transposed output): */

	initialize_fft_data(cin, local_N);
	fftwnd_mpi(plan_nd_mpi, 1, cin, FFTW_TRANSPOSED_ORDER);
	MPI_Barrier(MPI_COMM_WORLD);
	start_t = fftw_get_time();
	for (iter = 0; iter < max_iter; ++iter) {
	    initialize_fft_data(cin, local_N);
	    fftwnd_mpi(plan_nd_mpi, 1, cin, FFTW_TRANSPOSED_ORDER);
	}
	end_t = fftw_get_time();
	if (my_pe == 0) {
	printf("%17g", time2=fftw_time_to_sec(fftw_time_diff(fftw_time_diff(end_t,start_t),init_t)) *
                       time_scale);
	printf("%13g\n",time1/time2);
	fflush(stdout);
	}

	/* Done. */

	if (my_pe == 0)
	     fftwnd_destroy_plan(plan_nd);
	fftwnd_mpi_destroy_plan(plan_nd_mpi);

	if (my_pe != 0)
	     fftw_free(cin);
    }

    if (my_pe == 0)
	 fftw_free(cin);

    MPI_Finalize();

    return 0;
}

void initialize_fft_data(fftw_complex * arr, long n)
{
    long i;

    for (i = 0; i < n; i++) { /* initialize to some arbitrary values: */
	c_re(arr[i]) = 0.56923456;
	c_im(arr[i]) = 0.23858572;
    }
}
