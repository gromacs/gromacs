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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <fftw-int.h>
#include <rfftw_mpi.h>

#include <test_main.h>

#define my_printf if (io_okay) printf
#define my_fprintf if (io_okay) fprintf
#define my_fflush if (io_okay) fflush

int ncpus = 1;
int my_cpu = 0;

char fftw_prefix[] = "rfftw_mpi";

/*************************************************
 * Speed tests
 *************************************************/

#define MPI_TIME_FFT(fft,a,n,t) \
{ \
     double ts,te; \
     double total_t; \
     int iters = 1,iter; \
     zero_arr((n), (a)); \
     do { \
          MPI_Barrier(MPI_COMM_WORLD); \
          ts = MPI_Wtime(); \
          for (iter = 0; iter < iters; ++iter) fft; \
          te = MPI_Wtime(); \
          t = (total_t = (te - ts)) / iters; \
          iters *= 2; \
     } while (total_t < 2.0); \
}

void zero_arr(int n, fftw_real * a)
{
     int i;
     for (i = 0; i < n; ++i)
	  a[i] = 0.0;
}

void test_speed_aux(int n, fftw_direction dir, int flags, int specific)
{
     fftw_real *in, *out;
     fftw_plan plan;
     double t;
     fftw_time begin, end;

     return; /* one-dim transforms not supported yet in MPI */

     in = (fftw_real *) fftw_malloc(n * howmany_fields
				    * sizeof(fftw_real));
     out = (fftw_real *) fftw_malloc(n * howmany_fields
				     * sizeof(fftw_real));

     if (specific) {
	  begin = fftw_get_time();
	  plan = rfftw_create_plan_specific(n, dir,
					speed_flag | flags | wisdom_flag,
					    in, howmany_fields,
					    out, howmany_fields);
	  end = fftw_get_time();
     } else {
	  begin = fftw_get_time();
	  plan = rfftw_create_plan(n, dir, speed_flag | flags | wisdom_flag);
	  end = fftw_get_time();
     }
     CHECK(plan != NULL, "can't create plan");

     t = fftw_time_to_sec(fftw_time_diff(end, begin));
     WHEN_VERBOSE(2, printf("time for planner: %f s\n", t));

     WHEN_VERBOSE(2, rfftw_print_plan(plan));

     FFTW_TIME_FFT(rfftw(plan, howmany_fields,
			 in, howmany_fields, 1, out, howmany_fields, 1),
		   in, n * howmany_fields, t);

     rfftw_destroy_plan(plan);

     WHEN_VERBOSE(1, printf("time for one fft: %s", smart_sprint_time(t)));
     WHEN_VERBOSE(1, printf(" (%s/point)\n", smart_sprint_time(t / n)));
     WHEN_VERBOSE(1, printf("\"mflops\" = 5/2 (n log2 n) / (t in microseconds)"
			" = %f\n", 0.5 * howmany_fields * mflops(t, n)));

     fftw_free(in);
     fftw_free(out);

     WHEN_VERBOSE(1, printf("\n"));
}

void test_speed_nd_aux(struct size sz,
		       fftw_direction dir, int flags, int specific)
{
     int local_nx, local_x_start, local_ny_after_transpose,
          local_y_start_after_transpose, total_local_size;
     fftw_real *in, *work;
     rfftwnd_plan plan = 0;
     rfftwnd_mpi_plan mpi_plan;
     double t, t0 = 0.0;
     int i, N;

     if (sz.rank < 2)
          return;

     /* only bench in-place multi-dim transforms */
     flags |= FFTW_IN_PLACE;	

     N = 1;
     for (i = 0; i < sz.rank - 1; ++i)
	  N *= sz.narray[i];

     N *= (sz.narray[i] + 2);

     if (specific) {
	  return;
     } else {
          if (io_okay)
               plan = rfftwnd_create_plan(sz.rank, sz.narray,
                                        dir, speed_flag | flags | wisdom_flag);
          mpi_plan = rfftwnd_mpi_create_plan(MPI_COMM_WORLD, sz.rank,sz.narray,
                                        dir, speed_flag | flags | wisdom_flag);
     }
     CHECK(mpi_plan != NULL, "can't create plan");

     rfftwnd_mpi_local_sizes(mpi_plan, &local_nx, &local_x_start,
			     &local_ny_after_transpose,
			     &local_y_start_after_transpose,
			     &total_local_size);

     if (io_okay)
          in = (fftw_real *) fftw_malloc(N * howmany_fields *
					 sizeof(fftw_real));
     else
          in = (fftw_real *) fftw_malloc(total_local_size * howmany_fields *
					 sizeof(fftw_real));
     work = (fftw_real *) fftw_malloc(total_local_size * howmany_fields *
                                         sizeof(fftw_real));

     if (io_okay) {
	  if (dir == FFTW_REAL_TO_COMPLEX) {
	       FFTW_TIME_FFT(rfftwnd_real_to_complex(plan, howmany_fields,
						     in, howmany_fields, 1,
						     0, 0, 0),
			     in, N * howmany_fields, t0);
	  }
	  else {
	       FFTW_TIME_FFT(rfftwnd_complex_to_real(plan, howmany_fields,
						     (fftw_complex *) in,
						     howmany_fields, 1,
						     0, 0, 0),
			     in, N * howmany_fields, t0);
	  }
     }

     rfftwnd_destroy_plan(plan);

     WHEN_VERBOSE(1, my_printf("time for one fft (uniprocessor): %s\n",
                               smart_sprint_time(t0)));

     MPI_TIME_FFT(rfftwnd_mpi(mpi_plan, howmany_fields,
                             in, NULL, FFTW_NORMAL_ORDER),
                   in, total_local_size * howmany_fields, t);

     if (io_okay) {
          WHEN_VERBOSE(1, printf("NORMAL: time for one fft (%d cpus): %s",
                                 ncpus, smart_sprint_time(t)));
          WHEN_VERBOSE(1, printf(" (%s/point)\n", smart_sprint_time(t / N)));
          WHEN_VERBOSE(1, printf("NORMAL: \"mflops\" = 5/2 (N log2 N) / "
                                 "(t in microseconds)"
                                 " = %f\n", 0.5*howmany_fields*mflops(t, N)));
          WHEN_VERBOSE(1, printf("NORMAL: parallel speedup: %f\n", t0 / t));
     }

     MPI_TIME_FFT(rfftwnd_mpi(mpi_plan, howmany_fields,
                             in, NULL, FFTW_TRANSPOSED_ORDER),
                   in, total_local_size * howmany_fields, t);

     if (io_okay) {
          WHEN_VERBOSE(1, printf("TRANSP.: time for one fft (%d cpus): %s",
                                 ncpus, smart_sprint_time(t)));
          WHEN_VERBOSE(1, printf(" (%s/point)\n", smart_sprint_time(t / N)));
          WHEN_VERBOSE(1, printf("TRANSP.: \"mflops\" = 5/2 (N log2 N) / "
                                 "(t in microseconds)"
                                 " = %f\n", 0.5*howmany_fields*mflops(t, N)));
          WHEN_VERBOSE(1, printf("TRANSP.: parallel speedup: %f\n", t0 / t));
     }

     MPI_TIME_FFT(rfftwnd_mpi(mpi_plan, howmany_fields,
                             in, work, FFTW_NORMAL_ORDER),
                   in, total_local_size * howmany_fields, t);

     if (io_okay) {
          WHEN_VERBOSE(1, printf("NORMAL,w/WORK: time for one fft (%d cpus): %s",
                                 ncpus, smart_sprint_time(t)));
          WHEN_VERBOSE(1, printf(" (%s/point)\n", smart_sprint_time(t / N)));
          WHEN_VERBOSE(1, printf("NORMAL,w/WORK: \"mflops\" = 5/2 (N log2 N) / "
                                 "(t in microseconds)"
                                 " = %f\n", 0.5*howmany_fields*mflops(t, N)));
          WHEN_VERBOSE(1, printf("NORMAL,w/WORK: parallel speedup: %f\n",
				 t0 / t));
     }

     MPI_TIME_FFT(rfftwnd_mpi(mpi_plan, howmany_fields,
                             in, work, FFTW_TRANSPOSED_ORDER),
                   in, total_local_size * howmany_fields, t);

     if (io_okay) {
          WHEN_VERBOSE(1, printf("TRANSP.,w/WORK: time for one fft (%d cpus): %s",
                                 ncpus, smart_sprint_time(t)));
          WHEN_VERBOSE(1, printf(" (%s/point)\n", smart_sprint_time(t / N)));
          WHEN_VERBOSE(1, printf("TRANSP.,w/WORK: \"mflops\" = 5/2 (N log2 N) / "
                                 "(t in microseconds)"
                                 " = %f\n", 0.5*howmany_fields*mflops(t, N)));
          WHEN_VERBOSE(1, printf("TRANSP.,w/WORK: parallel speedup: %f\n",
				 t0 / t));
     }

     rfftwnd_mpi_destroy_plan(mpi_plan);

     fftw_free(in);
     fftw_free(work);

     WHEN_VERBOSE(1, my_printf("\n"));
}

/*************************************************
 * correctness tests
 *************************************************/

double compute_error(fftw_real * A, int astride,
                     fftw_real * B, int bstride, int n)
{
     /* compute the relative error */
     double error = 0.0;
     int i;

     for (i = 0; i < n; ++i) {
          double a;
          double mag;
          a = fabs(A[i * astride] - B[i * bstride]);
          mag = 0.5 * (fabs(A[i * astride]) + fabs(B[i * bstride]))+TOLERANCE;
          a /= mag;
          if (a > error)
               error = a;

#ifdef HAVE_ISNAN
          CHECK(!isnan(a), "NaN in answer");
#endif
     }
     return error;
}

void test_out_of_place(int n, int istride, int ostride,
		       int howmany, fftw_direction dir,
		       fftw_plan validated_plan, int specific)
{
     /* one-dim. out-of-place transforms will never be supported in MPI */
     WHEN_VERBOSE(2, my_printf("N/A\n"));
}

void test_in_place(int n, int istride,
		   int howmany, fftw_direction dir,
		   fftw_plan validated_plan, int specific)
{
     /* one-dim. transforms are not supported yet in MPI */
     WHEN_VERBOSE(2, my_printf("N/A\n"));
}

void test_out_of_place_both(int n, int istride, int ostride,
			    int howmany,
			    fftw_plan validated_plan_forward,
			    fftw_plan validated_plan_backward)
{
}

void test_in_place_both(int n, int istride, int howmany,
			fftw_plan validated_plan_forward,
			fftw_plan validated_plan_backward)
{
     WHEN_VERBOSE(2,
		  printf("TEST CORRECTNESS (in place, FFTW_FORWARD, %s) "
			 "n = %d  istride = %d  howmany = %d\n",
			 SPECIFICP(0),
			 n, istride, howmany));
     test_in_place(n, istride, howmany, FFTW_FORWARD,
		   validated_plan_forward, 0);
     
     WHEN_VERBOSE(2,
		  printf("TEST CORRECTNESS (in place, FFTW_BACKWARD, %s) "
			 "n = %d  istride = %d  howmany = %d\n",
			 SPECIFICP(0),
			 n, istride, howmany));
     test_in_place(n, istride, howmany, FFTW_BACKWARD,
		   validated_plan_backward, 0);
}

void test_correctness(int n)
{
}

/*************************************************
 * multi-dimensional correctness tests
 *************************************************/

void testnd_out_of_place(int rank, int *n, fftwnd_plan validated_plan)
{
}

void testnd_in_place(int rank, int *n, fftwnd_plan validated_plan,
		     int alternate_api, int specific)
{
     int local_nx, local_x_start, local_ny_after_transpose,
          local_y_start_after_transpose, total_local_size;
     int istride, ostride, howmany;
     int N, dim, i, j, k;
     int nc, nhc, nr;
     fftw_real *in1, *out3, *work = 0;
     fftw_complex *in2, *out1, *out2;
     rfftwnd_mpi_plan p = 0, ip = 0;
     int flags = measure_flag | wisdom_flag | FFTW_IN_PLACE;

     if (specific || rank < 2)
          return;

     if (coinflip())
	  flags |= FFTW_THREADSAFE;

     N = nc = nr = nhc = 1;
     for (dim = 0; dim < rank; ++dim)
	  N *= n[dim];
     if (rank > 0) {
	  nr = n[rank - 1];
	  nc = N / nr;
	  nhc = nr / 2 + 1;
     }

     if (alternate_api && (rank == 2 || rank == 3)) {
	  if (rank == 2) {
	       p = rfftw2d_mpi_create_plan(MPI_COMM_WORLD,
					   n[0], n[1], FFTW_REAL_TO_COMPLEX,
					   flags);
	       ip = rfftw2d_mpi_create_plan(MPI_COMM_WORLD,
					    n[0], n[1], FFTW_COMPLEX_TO_REAL,
					    flags);
	  } else {
	       p = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
					   n[0], n[1], n[2],
					   FFTW_REAL_TO_COMPLEX, flags);
	       ip = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
					    n[0], n[1], n[2],
					    FFTW_COMPLEX_TO_REAL, flags);
	  }
     }
     else {
	  p = rfftwnd_mpi_create_plan(MPI_COMM_WORLD,
				      rank, n, FFTW_REAL_TO_COMPLEX, flags);
	  ip = rfftwnd_mpi_create_plan(MPI_COMM_WORLD, 
				       rank, n, FFTW_COMPLEX_TO_REAL, flags);
     }

     CHECK(p != NULL && ip != NULL, "can't create plan");

     rfftwnd_mpi_local_sizes(p, &local_nx, &local_x_start,
			     &local_ny_after_transpose,
			     &local_y_start_after_transpose,
			     &total_local_size);
     
     in1 = (fftw_real *) fftw_malloc(total_local_size
				     * MAX_STRIDE * sizeof(fftw_real));
     if (coinflip()) {
          WHEN_VERBOSE(1, my_printf("w/work..."));
	  work = (fftw_real *) fftw_malloc(total_local_size 
					   * MAX_STRIDE * sizeof(fftw_real));
     }
     out3 = in1;
     out1 = (fftw_complex *) in1;
     in2 = (fftw_complex *) fftw_malloc(N * sizeof(fftw_complex));
     out2 = (fftw_complex *) fftw_malloc(N * sizeof(fftw_complex));
     for (i = 0; i < total_local_size * MAX_STRIDE; ++i)
	  out3[i] = 0;

     for (istride = 1; istride <= MAX_STRIDE; ++istride) {
	  /* generate random inputs */
	  for (i = 0; i < nc; ++i)
	       for (j = 0; j < nr; ++j) {
		    c_re(in2[i * nr + j]) = DRAND();
		    c_im(in2[i * nr + j]) = 0.0;
	       }

	  for (i = 0; i < local_nx * (nc / n[0]); ++i)
	       for (j = 0; j < nr; ++j) {
		    for (k = 0; k < istride; ++k)
			 in1[(i * nhc * 2 + j) * istride + k]
			     = c_re((in2 + local_x_start * (N/n[0]))
				     [i * nr + j]);
	       }


	  fftwnd(validated_plan, 1, in2, 1, 1, out2, 1, 1);

	  howmany = ostride = istride;

	  WHEN_VERBOSE(2, printf("\n    testing in-place stride %d...",
				 istride));

	  rfftwnd_mpi(p, howmany, in1, work, FFTW_NORMAL_ORDER);

	  for (i = 0; i < local_nx * (nc / n[0]); ++i)
	       for (k = 0; k < howmany; ++k)
		    CHECK(compute_error_complex(out1 + i * nhc * ostride + k,
						ostride,
						out2 + local_x_start*(N/n[0])
						+ i * nr, 1,
						nhc) < TOLERANCE,
			  "in-place (r2c): wrong answer");

	  rfftwnd_mpi(ip, howmany, in1, work, FFTW_NORMAL_ORDER);

	  for (i = 0; i < total_local_size * istride; ++i)
	       out3[i] *= 1.0 / N;

	  for (i = 0; i < local_nx * (nc / n[0]); ++i)
	       for (k = 0; k < howmany; ++k)
		    CHECK(compute_error(out3 + i * nhc * 2 * istride + k,
					istride,
					(fftw_real *)
					(in2 + local_x_start*(N/n[0])
					 + i * nr), 2,
					nr) < TOLERANCE,
			  "in-place (c2r): wrong answer (check 2)");
     }

     rfftwnd_mpi_destroy_plan(p);
     rfftwnd_mpi_destroy_plan(ip);

     fftw_free(work);
     fftw_free(out2);
     fftw_free(in2);
     fftw_free(in1);
}

void testnd_correctness(struct size sz, fftw_direction dir,
			int alt_api, int specific, int force_buf)
{
     fftwnd_plan validated_plan;

     if (dir != FFTW_FORWARD)
	  return;
     if (force_buf)
	  return;

     validated_plan = fftwnd_create_plan(sz.rank, sz.narray, 
					 dir, measure_flag | wisdom_flag);
     CHECK(validated_plan != NULL, "can't create plan");

     testnd_in_place(sz.rank, sz.narray,
		     validated_plan, alt_api, specific);

     fftwnd_destroy_plan(validated_plan);
}

/*************************************************
 * planner tests
 *************************************************/

void test_planner(int rank)
{
     /*
      * create and destroy many plans, at random.  Check the
      * garbage-collecting allocator of twiddle factors
      */
     int i, dim;
     int r, s;
     rfftwnd_mpi_plan pnd[PLANNER_TEST_SIZE];
     int *narr, maxdim;

     chk_mem_leak = 0;
     verbose--;

     please_wait();
     if (rank < 2)
          rank = 2;

     narr = (int *) fftw_malloc(rank * sizeof(int));

     for (i = 0; i < PLANNER_TEST_SIZE; ++i) {
          pnd[i] = (rfftwnd_mpi_plan) 0;
     }

     maxdim = (int) pow(8192, 1.0/rank);

     for (i = 0; i < PLANNER_TEST_SIZE * PLANNER_TEST_SIZE; ++i) {
          r = rand();
          if (r < 0)
               r = -r;
          r = r % PLANNER_TEST_SIZE;

          for (dim = 0; dim < rank; ++dim) {
               do {
                    s = rand();
                    if (s < 0)
                         s = -s;
                    s = s % maxdim + 1;
               } while (s == 0);
               narr[dim] = s;
          }

	  if (rank > 1) {
	       if (pnd[r])
		    rfftwnd_mpi_destroy_plan(pnd[r]);
	       
	       pnd[r] = rfftwnd_mpi_create_plan(MPI_COMM_WORLD, rank, narr,
					       random_dir(), measure_flag |
					       wisdom_flag);
	  }

          if (i % (PLANNER_TEST_SIZE * PLANNER_TEST_SIZE / 20) == 0) {
               WHEN_VERBOSE(0, my_printf("test planner: so far so good\n"));
               WHEN_VERBOSE(0, my_printf("test planner: iteration %d"
					 " out of %d\n",
                              i, PLANNER_TEST_SIZE * PLANNER_TEST_SIZE));
          }
     }

     for (i = 0; i < PLANNER_TEST_SIZE; ++i) {
          if (pnd[i])
               rfftwnd_mpi_destroy_plan(pnd[i]);
     }

     fftw_free(narr);
     verbose++;
     chk_mem_leak = 1;
}

/*************************************************
 * test initialization
 *************************************************/

void test_init(int *argc, char **argv)
{
     unsigned int seed;

     MPI_Init(argc,&argv);
     MPI_Comm_size(MPI_COMM_WORLD,&ncpus);
     MPI_Comm_rank(MPI_COMM_WORLD,&my_cpu);

     /* Only process 0 gets to do I/O: */
     io_okay = my_cpu == 0;

     /* Make sure all processes use the same seed for random numbers: */
     seed = time(NULL);
     MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
     srand(seed);

     fftw_die_hook = fftw_mpi_die; /* call MPI_Abort on failure */
}

void test_finish(void)
{
     MPI_Finalize();
}

void enter_paranoid_mode(void)
{
}

/* in MPI, only process 0 is guaranteed to have access to the argument list */
int get_option(int argc, char **argv, char *argval, int argval_maxlen)
{
     int c;
     int arglen;

     if (io_okay) {
	  c = default_get_option(argc,argv,argval,argval_maxlen);
	  arglen = strlen(argval) + 1;
     }

     MPI_Bcast(&c, 1, MPI_INT, 0, MPI_COMM_WORLD);
     MPI_Bcast(&arglen, 1, MPI_INT, 0, MPI_COMM_WORLD);
     MPI_Bcast(argval, arglen, MPI_CHAR, 0, MPI_COMM_WORLD);

     return c;
}
