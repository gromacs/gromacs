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

/*
 * fftw_test.c : test program for complex-complex transforms
 */

/* $Id$ */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <fftw-int.h>
#include <fftw_threads.h>

#include <test_main.h>

char fftw_prefix[] = "fftw_threads";
int nthreads = 1;

/*************************************************
 * Speed tests
 *************************************************/

void zero_arr(int n, fftw_complex * a)
{
     int i;
     for (i = 0; i < n; ++i)
	  c_re(a[i]) = c_im(a[i]) = 0.0;
}

void test_speed_aux(int n, fftw_direction dir, int flags, int specific)
{
     fftw_complex *in, *out;
     fftw_plan plan;
     double t, t0;
     fftw_time begin, end;

     in = (fftw_complex *) fftw_malloc(n * howmany_fields
				       * sizeof(fftw_complex));
     out = (fftw_complex *) fftw_malloc(n * howmany_fields
					* sizeof(fftw_complex));

     if (specific) {
	  begin = fftw_get_time();
	  plan = fftw_create_plan_specific(n, dir,
					speed_flag | flags | wisdom_flag,
					   in, howmany_fields,
					   out, howmany_fields);
	  end = fftw_get_time();
     } else {
	  begin = fftw_get_time();
	  plan = fftw_create_plan(n, dir, speed_flag | flags | wisdom_flag);
	  end = fftw_get_time();
     }
     CHECK(plan != NULL, "can't create plan");

     t = fftw_time_to_sec(fftw_time_diff(end, begin));
     WHEN_VERBOSE(2, printf("time for planner: %f s\n", t));

     WHEN_VERBOSE(2, fftw_print_plan(plan));

     FFTW_TIME_FFT(fftw(plan, howmany_fields,
			in, howmany_fields, 1, out, howmany_fields, 1),
		   in, n * howmany_fields, t0);
     
     FFTW_TIME_FFT(fftw_threads(nthreads, plan, howmany_fields,
				in, howmany_fields, 1, out, howmany_fields, 1),
		   in, n * howmany_fields, t);

     fftw_destroy_plan(plan);

     WHEN_VERBOSE(1, printf("time for one fft (uniprocessor): %s\n", smart_sprint_time(t0)));
     WHEN_VERBOSE(1, printf("time for one fft (%d threads): %s", nthreads, smart_sprint_time(t)));
     WHEN_VERBOSE(1, printf(" (%s/point)\n", smart_sprint_time(t / n)));
     WHEN_VERBOSE(1, printf("\"mflops\" = 5 (n log2 n) / (t in microseconds)"
			    " = %f\n", howmany_fields * mflops(t, n)));
     WHEN_VERBOSE(1, printf("parallel speedup: %f\n", t0 / t));

     fftw_free(in);
     fftw_free(out);

     WHEN_VERBOSE(1, printf("\n"));
}

void test_speed_nd_aux(struct size sz,
		       fftw_direction dir, int flags, int specific)
{
     fftw_complex *in;
     fftwnd_plan plan;
     double t, t0;
     fftw_time begin, end;
     int i, N;

     /* only bench in-place multi-dim transforms */
     flags |= FFTW_IN_PLACE;	

     N = 1;
     for (i = 0; i < sz.rank; ++i)
	  N *= (sz.narray[i]);

     in = (fftw_complex *) fftw_malloc(N * howmany_fields *
				       sizeof(fftw_complex));

     if (specific) {
	  begin = fftw_get_time();
	  plan = fftwnd_create_plan_specific(sz.rank, sz.narray, dir,
					speed_flag | flags | wisdom_flag,
					     in, howmany_fields, 0, 1);
     } else {
	  begin = fftw_get_time();
	  plan = fftwnd_create_plan(sz.rank, sz.narray,
				  dir, speed_flag | flags | wisdom_flag);
     }
     end = fftw_get_time();
     CHECK(plan != NULL, "can't create plan");

     t = fftw_time_to_sec(fftw_time_diff(end, begin));
     WHEN_VERBOSE(2, printf("time for planner: %f s\n", t));

     WHEN_VERBOSE(2, printf("\n"));
     WHEN_VERBOSE(2, (fftwnd_print_plan(plan)));
     WHEN_VERBOSE(2, printf("\n"));

     FFTW_TIME_FFT(fftwnd(plan, howmany_fields,
			  in, howmany_fields, 1, 0, 0, 0),
		   in, N * howmany_fields, t0);

     FFTW_TIME_FFT(fftwnd_threads(nthreads, plan, howmany_fields,
				  in, howmany_fields, 1, 0, 0, 0),
		   in, N * howmany_fields, t);

     fftwnd_destroy_plan(plan);

     WHEN_VERBOSE(1, printf("time for one fft (uniprocessor): %s\n", smart_sprint_time(t0)));
     WHEN_VERBOSE(1, printf("time for one fft (%d threads): %s", nthreads, smart_sprint_time(t)));
     WHEN_VERBOSE(1, printf(" (%s/point)\n", smart_sprint_time(t / N)));
     WHEN_VERBOSE(1, printf("\"mflops\" = 5 (N log2 N) / (t in microseconds)"
			    " = %f\n", howmany_fields * mflops(t, N)));
     WHEN_VERBOSE(1, printf("parallel speedup: %f\n", t0 / t));

     fftw_free(in);

     WHEN_VERBOSE(1, printf("\n"));
}

/*************************************************
 * correctness tests
 *************************************************/

void test_out_of_place(int n, int istride, int ostride,
		       int howmany, fftw_direction dir,
		       fftw_plan validated_plan,
		       int specific)
{
     fftw_complex *in1, *in2, *out1, *out2;
     fftw_plan plan;
     int i, j;
     int flags = measure_flag | wisdom_flag | FFTW_OUT_OF_PLACE;

     if (coinflip())
	  flags |= FFTW_THREADSAFE;

     in1 = (fftw_complex *) fftw_malloc(istride * n * sizeof(fftw_complex) * howmany);
     in2 = (fftw_complex *) fftw_malloc(n * sizeof(fftw_complex) * howmany);
     out1 = (fftw_complex *) fftw_malloc(ostride * n * sizeof(fftw_complex) * howmany);
     out2 = (fftw_complex *) fftw_malloc(n * sizeof(fftw_complex) * howmany);

     if (!specific)
	  plan = fftw_create_plan(n, dir, flags);
     else
	  plan = fftw_create_plan_specific(n, dir,
					   flags,
					   in1, istride, out1, ostride);

     /* generate random inputs */
     for (i = 0; i < n * howmany; ++i) {
	  c_re(in1[i * istride]) = c_re(in2[i]) = DRAND();
	  c_im(in1[i * istride]) = c_im(in2[i]) = DRAND();
     }

     /* 
      * fill in other positions of the array, to make sure that
      * fftw doesn't overwrite them 
      */
     for (j = 1; j < istride; ++j)
	  for (i = 0; i < n * howmany; ++i) {
	       c_re(in1[i * istride + j]) = i * istride + j;
	       c_im(in1[i * istride + j]) = i * istride - j;
	  }

     for (j = 1; j < ostride; ++j)
	  for (i = 0; i < n * howmany; ++i) {
	       c_re(out1[i * ostride + j]) = -i * ostride + j;
	       c_im(out1[i * ostride + j]) = -i * ostride - j;
	  }

     CHECK(plan != NULL, "can't create plan");
     WHEN_VERBOSE(2, fftw_print_plan(plan));

     /* fft-ize */
     if (howmany != 1 || istride != 1 || ostride != 1 || coinflip())
	  fftw_threads(nthreads, plan, howmany, in1, istride, n * istride,
		       out1, ostride, n * ostride);
     else
	  fftw_threads_one(nthreads, plan, in1, out1);

     fftw_destroy_plan(plan);

     /* check for overwriting */
     for (j = 1; j < istride; ++j)
	  for (i = 0; i < n * howmany; ++i)
	       CHECK(c_re(in1[i * istride + j]) == i * istride + j &&
		     c_im(in1[i * istride + j]) == i * istride - j,
		     "input has been overwritten");
     for (j = 1; j < ostride; ++j)
	  for (i = 0; i < n * howmany; ++i)
	       CHECK(c_re(out1[i * ostride + j]) == -i * ostride + j &&
		     c_im(out1[i * ostride + j]) == -i * ostride - j,
		     "output has been overwritten");

     for (i = 0; i < howmany; ++i) {
	  fftw(validated_plan, 1, in2 + n * i, 1, n, out2 + n * i, 1, n);
     }

     CHECK(compute_error_complex(out1, ostride, out2, 1, n * howmany) < TOLERANCE,
	   "test_out_of_place: wrong answer");
     WHEN_VERBOSE(2, printf("OK\n"));

     fftw_free(in1);
     fftw_free(in2);
     fftw_free(out1);
     fftw_free(out2);
}

void test_in_place(int n, int istride, int howmany, fftw_direction dir,
		   fftw_plan validated_plan, int specific)
{
     fftw_complex *in1, *in2, *out2;
     fftw_plan plan;
     int i, j;
     int flags = measure_flag | wisdom_flag | FFTW_IN_PLACE;

     if (coinflip())
	  flags |= FFTW_THREADSAFE;

     in1 = (fftw_complex *) fftw_malloc(istride * n * sizeof(fftw_complex) * howmany);
     in2 = (fftw_complex *) fftw_malloc(n * sizeof(fftw_complex) * howmany);
     out2 = (fftw_complex *) fftw_malloc(n * sizeof(fftw_complex) * howmany);

     if (!specific)
	  plan = fftw_create_plan(n, dir, flags);
     else
	  plan = fftw_create_plan_specific(n, dir, flags,
					   in1, istride,
					   (fftw_complex *) NULL, 0);

     /* generate random inputs */
     for (i = 0; i < n * howmany; ++i) {
	  c_re(in1[i * istride]) = c_re(in2[i]) = DRAND();
	  c_im(in1[i * istride]) = c_im(in2[i]) = DRAND();
     }

     /* 
      * fill in other positions of the array, to make sure that
      * fftw doesn't overwrite them 
      */
     for (j = 1; j < istride; ++j)
	  for (i = 0; i < n * howmany; ++i) {
	       c_re(in1[i * istride + j]) = i * istride + j;
	       c_im(in1[i * istride + j]) = i * istride - j;
	  }
     CHECK(plan != NULL, "can't create plan");
     WHEN_VERBOSE(2, fftw_print_plan(plan));

     /* fft-ize */
     if (howmany != 1 || istride != 1 || coinflip())
	  fftw_threads(nthreads, plan, howmany, in1, istride, n * istride,
		       (fftw_complex *) NULL, 0, 0);
     else
	  fftw_threads_one(nthreads, plan, in1, NULL);

     fftw_destroy_plan(plan);

     /* check for overwriting */
     for (j = 1; j < istride; ++j)
	  for (i = 0; i < n * howmany; ++i)
	       CHECK(c_re(in1[i * istride + j]) == i * istride + j &&
		     c_im(in1[i * istride + j]) == i * istride - j,
		     "input has been overwritten");

     for (i = 0; i < howmany; ++i) {
	  fftw(validated_plan, 1, in2 + n * i, 1, n, out2 + n * i, 1, n);
     }

     CHECK(compute_error_complex(in1, istride, out2, 1, n * howmany) < TOLERANCE,
	   "test_in_place: wrong answer");
     WHEN_VERBOSE(2, printf("OK\n"));

     fftw_free(in1);
     fftw_free(in2);
     fftw_free(out2);
}

void test_out_of_place_both(int n, int istride, int ostride,
			    int howmany,
			    fftw_plan validated_plan_forward,
			    fftw_plan validated_plan_backward)
{
     int specific;

     for (specific = 0; specific <= 1; ++specific) {
	  WHEN_VERBOSE(2,
	       printf("TEST CORRECTNESS (out of place, FFTW_FORWARD, %s)"
		   " n = %d  istride = %d  ostride = %d  howmany = %d\n",
		      SPECIFICP(specific),
		      n, istride, ostride, howmany));
	  test_out_of_place(n, istride, ostride, howmany, FFTW_FORWARD,
			    validated_plan_forward, specific);

	  WHEN_VERBOSE(2,
	      printf("TEST CORRECTNESS (out of place, FFTW_BACKWARD, %s)"
		   " n = %d  istride = %d  ostride = %d  howmany = %d\n",
		     SPECIFICP(specific),
		     n, istride, ostride, howmany));
	  test_out_of_place(n, istride, ostride, howmany, FFTW_BACKWARD,
			    validated_plan_backward, specific);
     }
}

void test_in_place_both(int n, int istride, int howmany,
			fftw_plan validated_plan_forward,
			fftw_plan validated_plan_backward)
{
     int specific;

     for (specific = 0; specific <= 1; ++specific) {
	  WHEN_VERBOSE(2,
		  printf("TEST CORRECTNESS (in place, FFTW_FORWARD, %s) "
			 "n = %d  istride = %d  howmany = %d\n",
			 SPECIFICP(specific),
			 n, istride, howmany));
	  test_in_place(n, istride, howmany, FFTW_FORWARD,
			validated_plan_forward, specific);

	  WHEN_VERBOSE(2,
		 printf("TEST CORRECTNESS (in place, FFTW_BACKWARD, %s) "
			"n = %d  istride = %d  howmany = %d\n",
			SPECIFICP(specific),
			n, istride, howmany));
	  test_in_place(n, istride, howmany, FFTW_BACKWARD,
			validated_plan_backward, specific);
     }
}

void test_correctness(int n)
{
     int istride, ostride, howmany;
     fftw_plan validated_plan_forward, validated_plan_backward;

     WHEN_VERBOSE(1,
		  printf("Testing correctness for n = %d...", n);
		  fflush(stdout));

     /* produce a good plan */
     validated_plan_forward =
	 fftw_create_plan(n, FFTW_FORWARD, measure_flag | wisdom_flag);
     validated_plan_backward =
	 fftw_create_plan(n, FFTW_BACKWARD, measure_flag | wisdom_flag);

     for (istride = 1; istride <= MAX_STRIDE; ++istride)
	  for (ostride = 1; ostride <= MAX_STRIDE; ++ostride)
	       for (howmany = 1; howmany <= MAX_HOWMANY; ++howmany)
		    test_out_of_place_both(n, istride, ostride, howmany,
					   validated_plan_forward,
					   validated_plan_backward);

     for (istride = 1; istride <= MAX_STRIDE; ++istride)
	  for (howmany = 1; howmany <= MAX_HOWMANY; ++howmany)
	       test_in_place_both(n, istride, howmany,
				  validated_plan_forward,
				  validated_plan_backward);

     fftw_destroy_plan(validated_plan_forward);
     fftw_destroy_plan(validated_plan_backward);

     if (!(wisdom_flag & FFTW_USE_WISDOM) && chk_mem_leak)
	  fftw_check_memory_leaks();

     WHEN_VERBOSE(1, printf("OK\n"));
}

/*************************************************
 * multi-dimensional correctness tests
 *************************************************/

void testnd_out_of_place(int rank, int *n, fftw_direction dir,
			 fftwnd_plan validated_plan)
{
     int istride, ostride;
     int N, dim, i;
     fftw_complex *in1, *in2, *out1, *out2;
     fftwnd_plan p;
     int flags = measure_flag | wisdom_flag;

     if (coinflip())
	  flags |= FFTW_THREADSAFE;

     N = 1;
     for (dim = 0; dim < rank; ++dim)
	  N *= n[dim];

     in1 = (fftw_complex *) fftw_malloc(N * MAX_STRIDE * sizeof(fftw_complex));
     out1 = (fftw_complex *) fftw_malloc(N * MAX_STRIDE * sizeof(fftw_complex));
     in2 = (fftw_complex *) fftw_malloc(N * sizeof(fftw_complex));
     out2 = (fftw_complex *) fftw_malloc(N * sizeof(fftw_complex));

     p = fftwnd_create_plan(rank, n, dir, flags);

     for (istride = 1; istride <= MAX_STRIDE; ++istride) {
	  /* generate random inputs */
	  for (i = 0; i < N; ++i) {
	       int j;
	       c_re(in2[i]) = DRAND();
	       c_im(in2[i]) = DRAND();
	       for (j = 0; j < istride; ++j) {
		    c_re(in1[i * istride + j]) = c_re(in2[i]);
		    c_im(in1[i * istride + j]) = c_im(in2[i]);
	       }
	  }

	  for (ostride = 1; ostride <= MAX_STRIDE; ++ostride) {
	       int howmany = (istride < ostride) ? istride : ostride;

	       if (howmany != 1 || istride != 1 || ostride != 1 || coinflip())
		    fftwnd_threads(nthreads, p, howmany, in1, istride, 1,
				   out1, ostride, 1);
	       else
		    fftwnd_threads_one(nthreads, p, in1, out1);

	       fftwnd(validated_plan, 1, in2, 1, 1, out2, 1, 1);

	       for (i = 0; i < howmany; ++i)
		    CHECK(compute_error_complex(out1 + i, ostride, out2, 1, N)
			  < TOLERANCE,
			  "testnd_out_of_place: wrong answer");
	  }
     }

     fftwnd_destroy_plan(p);

     fftw_free(out2);
     fftw_free(in2);
     fftw_free(out1);
     fftw_free(in1);
}

void testnd_in_place(int rank, int *n, fftw_direction dir,
		     fftwnd_plan validated_plan,
		     int alternate_api, int specific, int force_buffered)
{
     int istride;
     int N, dim, i;
     fftw_complex *in1, *in2, *out2;
     fftwnd_plan p;
     int flags = measure_flag | wisdom_flag | FFTW_IN_PLACE;

     if (coinflip())
	  flags |= FFTW_THREADSAFE;

     if (force_buffered)
	  flags |= FFTWND_FORCE_BUFFERED;

     N = 1;
     for (dim = 0; dim < rank; ++dim)
	  N *= n[dim];

     in1 = (fftw_complex *) fftw_malloc(N * MAX_STRIDE * sizeof(fftw_complex));
     in2 = (fftw_complex *) fftw_malloc(N * sizeof(fftw_complex));
     out2 = (fftw_complex *) fftw_malloc(N * sizeof(fftw_complex));

     if (!specific) {
	  if (alternate_api && (rank == 2 || rank == 3)) {
	       if (rank == 2)
		    p = fftw2d_create_plan(n[0], n[1], dir, flags);
	       else
		    p = fftw3d_create_plan(n[0], n[1], n[2], dir, flags);
	  } else		/* standard api */
	       p = fftwnd_create_plan(rank, n, dir, flags);
     } else {			/* specific plan creation */
	  if (alternate_api && (rank == 2 || rank == 3)) {
	       if (rank == 2)
		    p = fftw2d_create_plan_specific(n[0], n[1], dir, flags,
						    in1, 1,
					       (fftw_complex *) NULL, 1);
	       else
		    p = fftw3d_create_plan_specific(n[0], n[1], n[2], dir, flags,
						    in1, 1,
					       (fftw_complex *) NULL, 1);
	  } else		/* standard api */
	       p = fftwnd_create_plan_specific(rank, n, dir, flags,
					       in1, 1,
					       (fftw_complex *) NULL, 1);

     }

     for (istride = 1; istride <= MAX_STRIDE; ++istride) {
	  /* 
	   * generate random inputs */
	  for (i = 0; i < N; ++i) {
	       int j;
	       c_re(in2[i]) = DRAND();
	       c_im(in2[i]) = DRAND();
	       for (j = 0; j < istride; ++j) {
		    c_re(in1[i * istride + j]) = c_re(in2[i]);
		    c_im(in1[i * istride + j]) = c_im(in2[i]);
	       }
	  }

	  if (istride != 1 || istride != 1 || coinflip())
	       fftwnd_threads(nthreads, p, istride, in1, istride, 1,
			      (fftw_complex *) NULL, 1, 1);
	  else
	       fftwnd_threads_one(nthreads, p, in1, NULL);

	  fftwnd(validated_plan, 1, in2, 1, 1, out2, 1, 1);

	  for (i = 0; i < istride; ++i)
	       CHECK(compute_error_complex(in1 + i, istride, out2, 1, N) < TOLERANCE,
		     "testnd_in_place: wrong answer");
     }

     fftwnd_destroy_plan(p);

     fftw_free(out2);
     fftw_free(in2);
     fftw_free(in1);
}

void testnd_correctness(struct size sz, fftw_direction dir,
			int alt_api, int specific, int force_buf)
{
     fftwnd_plan validated_plan;

     validated_plan = fftwnd_create_plan(sz.rank, sz.narray,
					 dir, measure_flag | wisdom_flag);

     testnd_out_of_place(sz.rank, sz.narray, dir, validated_plan);
     testnd_in_place(sz.rank, sz.narray, dir, validated_plan, alt_api, specific, force_buf);

     fftwnd_destroy_plan(validated_plan);
}

/*************************************************
 * planner tests
 *************************************************/

void test_planner(int rank)
{
     WHEN_VERBOSE(1, printf("Use fftw_test to test the planner.\n"););
}

/*************************************************
 * test initialization
 *************************************************/

void test_init(int *argc, char **argv)
{
     int i;

     if (*argc >= 2)
	  nthreads = atoi(argv[1]);

     if (nthreads <= 0) {
	  fprintf(stderr, "Usage: fftw_threads_test <nthreads> [ options ]\n");
	  exit(EXIT_FAILURE);
     }
     for (i = 2; i < *argc; ++i)
	  argv[i - 1] = argv[i];
     *argc -= 1;

     if (fftw_threads_init()) {
	  fprintf(stderr, "Error initializing threads!");
	  exit(EXIT_FAILURE);
     }
}

void test_finish(void)
{
}

void enter_paranoid_mode(void)
{
}

int get_option(int argc, char **argv, char *argval, int argval_maxlen)
{
     return default_get_option(argc,argv,argval,argval_maxlen);
}
