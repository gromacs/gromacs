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

#include <fftw_threads-int.h>
#include <rfftw_threads.h>

#include <test_main.h>

char fftw_prefix[] = "rfftw_threads";
int nthreads = 1;

/*************************************************
 * Speed tests
 *************************************************/

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
     double t, t0;
     fftw_time begin, end;

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
		   in, n * howmany_fields, t0);

     FFTW_TIME_FFT(rfftw_threads(nthreads,plan, howmany_fields,
			 in, howmany_fields, 1, out, howmany_fields, 1),
		   in, n * howmany_fields, t);

     rfftw_destroy_plan(plan);

     WHEN_VERBOSE(1, printf("time for one fft (uniprocessor): %s\n", smart_sprint_time(t0)));
     WHEN_VERBOSE(1, printf("time for one fft (%d threads): %s", nthreads, smart_sprint_time(t)));
     WHEN_VERBOSE(1, printf(" (%s/point)\n", smart_sprint_time(t / n)));
     WHEN_VERBOSE(1, printf("\"mflops\" = 5/2 (n log2 n) / (t in microseconds)"
			" = %f\n", 0.5 * howmany_fields * mflops(t, n)));
     WHEN_VERBOSE(1, printf("parallel speedup: %f\n", t0 / t));

     fftw_free(in);
     fftw_free(out);

     WHEN_VERBOSE(1, printf("\n"));
}

void test_speed_nd_aux(struct size sz,
		       fftw_direction dir, int flags, int specific)
{
     fftw_real *in;
     fftwnd_plan plan;
     double t, t0;
     fftw_time begin, end;
     int i, N;

     /* only bench in-place multi-dim transforms */
     flags |= FFTW_IN_PLACE;	

     N = 1;
     for (i = 0; i < sz.rank - 1; ++i)
	  N *= sz.narray[i];

     N *= (sz.narray[i] + 2);

     in = (fftw_real *) fftw_malloc(N * howmany_fields * sizeof(fftw_real));

     if (specific) {
	  begin = fftw_get_time();
	  plan = rfftwnd_create_plan_specific(sz.rank, sz.narray, dir,
					speed_flag | flags | wisdom_flag,
					      in, howmany_fields, 0, 1);
     } else {
	  begin = fftw_get_time();
	  plan = rfftwnd_create_plan(sz.rank, sz.narray,
				  dir, speed_flag | flags | wisdom_flag);
     }
     end = fftw_get_time();
     CHECK(plan != NULL, "can't create plan");

     t = fftw_time_to_sec(fftw_time_diff(end, begin));
     WHEN_VERBOSE(2, printf("time for planner: %f s\n", t));

     WHEN_VERBOSE(2, printf("\n"));
     WHEN_VERBOSE(2, (rfftwnd_print_plan(plan)));
     WHEN_VERBOSE(2, printf("\n"));

     if (dir == FFTW_REAL_TO_COMPLEX) {
	  FFTW_TIME_FFT(rfftwnd_real_to_complex(plan, howmany_fields,
						in, howmany_fields, 1,
						0, 0, 0),
			in, N * howmany_fields, t0);
	  FFTW_TIME_FFT(rfftwnd_threads_real_to_complex(nthreads, plan, howmany_fields,
						in, howmany_fields, 1,
						0, 0, 0),
			in, N * howmany_fields, t);
     } else {
	  FFTW_TIME_FFT(rfftwnd_complex_to_real(plan, howmany_fields,
						(fftw_complex *) in,
						howmany_fields, 1,
						0, 0, 0),
			in, N * howmany_fields, t0);
	  FFTW_TIME_FFT(rfftwnd_threads_complex_to_real(nthreads, plan, howmany_fields,
						(fftw_complex *) in,
						howmany_fields, 1,
						0, 0, 0),
			in, N * howmany_fields, t);
     }

     rfftwnd_destroy_plan(plan);

     WHEN_VERBOSE(1, printf("time for one fft (uniprocessor): %s\n", smart_sprint_time(t0)));
     WHEN_VERBOSE(1, printf("time for one fft (%d threads): %s", nthreads, smart_sprint_time(t)));
     WHEN_VERBOSE(1, printf(" (%s/point)\n", smart_sprint_time(t / N)));
     WHEN_VERBOSE(1, printf("\"mflops\" = 5/2 (N log2 N) / (t in microseconds)"
			" = %f\n", 0.5 * howmany_fields * mflops(t, N)));
     WHEN_VERBOSE(1, printf("parallel speedup: %f\n", t0 / t));

     fftw_free(in);

     WHEN_VERBOSE(1, printf("\n"));
}

/*************************************************
 * correctness tests
 *************************************************/

void fill_random(fftw_real * a, int n, int stride)
{
     int i;

     /* generate random inputs */
     for (i = 0; i < n; ++i)
	  a[i * stride] = DRAND();
}

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
	  mag = 0.5 * (fabs(A[i * astride]) + fabs(B[i * bstride])) + TOLERANCE;

	  a /= mag;
	  if (a > error)
	       error = a;

#ifdef HAVE_ISNAN
	  CHECK(!isnan(a), "NaN in answer");
#endif
     }
     return error;
}

void array_compare(fftw_real * A, fftw_real * B, int n)
{
     CHECK(compute_error(A, 1, B, 1, n) < TOLERANCE,
	   "failure in RFFTW verification");
}

void test_out_of_place(int n, int istride, int ostride,
		       int howmany, fftw_direction dir,
		       fftw_plan validated_plan, int specific)
{
     fftw_complex *in2, *out2;
     fftw_real *in1, *out1, *out3;
     fftw_plan plan;
     int i, j;
     int flags = measure_flag | wisdom_flag;

     if (coinflip())
	  flags |= FFTW_THREADSAFE;

     in1 = (fftw_real *) fftw_malloc(istride * n * sizeof(fftw_real) * howmany);
     in2 = (fftw_complex *) fftw_malloc(n * sizeof(fftw_complex));
     out1 = (fftw_real *) fftw_malloc(ostride * n * sizeof(fftw_real) * howmany);
     out2 = (fftw_complex *) fftw_malloc(n * sizeof(fftw_complex));
     out3 = (fftw_real *) fftw_malloc(n * sizeof(fftw_real));

     if (!specific)
	  plan = rfftw_create_plan(n, dir, flags);
     else
	  plan = rfftw_create_plan_specific(n, dir, flags,
					    in1, istride, out1, ostride);
     CHECK(plan != NULL, "can't create plan");

     /* generate random inputs */
     fill_random(in1, n, istride);
     for (j = 1; j < howmany; ++j)
	  for (i = 0; i < n; ++i)
	       in1[(j * n + i) * istride] = in1[i * istride];

     /* copy random inputs to complex array for comparison with fftw: */
     if (dir == FFTW_REAL_TO_COMPLEX)
	  for (i = 0; i < n; ++i) {
	       c_re(in2[i]) = in1[i * istride];
	       c_im(in2[i]) = 0.0;
     } else {
	  int n2 = (n + 1) / 2;
	  c_re(in2[0]) = in1[0];
	  c_im(in2[0]) = 0.0;
	  for (i = 1; i < n2; ++i) {
	       c_re(in2[i]) = in1[i * istride];
	       c_im(in2[i]) = in1[(n - i) * istride];
	  }
	  if (n2 * 2 == n) {
	       c_re(in2[n2]) = in1[n2 * istride];
	       c_im(in2[n2]) = 0.0;
	       ++i;
	  }
	  for (; i < n; ++i) {
	       c_re(in2[i]) = c_re(in2[n - i]);
	       c_im(in2[i]) = -c_im(in2[n - i]);
	  }
     }

     /* 
      * fill in other positions of the array, to make sure that
      * rfftw doesn't overwrite them 
      */
     for (j = 1; j < istride; ++j)
	  for (i = 0; i < n * howmany; ++i)
	       in1[i * istride + j] = i * istride + j;

     for (j = 1; j < ostride; ++j)
	  for (i = 0; i < n * howmany; ++i)
	       out1[i * ostride + j] = -i * ostride + j;

     WHEN_VERBOSE(2, rfftw_print_plan(plan));

     /* fft-ize */
     if (howmany != 1 || istride != 1 || ostride != 1 || coinflip())
	  rfftw_threads(nthreads, plan, howmany, in1, istride, n * istride,
			out1, ostride, n * ostride);
     else
	  rfftw_threads_one(nthreads, plan, in1, out1);

     rfftw_destroy_plan(plan);

     /* check for overwriting */
     for (j = 1; j < istride; ++j)
	  for (i = 0; i < n * howmany; ++i)
	       CHECK(in1[i * istride + j] == i * istride + j,
		     "input has been overwritten");
     for (j = 1; j < ostride; ++j)
	  for (i = 0; i < n * howmany; ++i)
	       CHECK(out1[i * ostride + j] == -i * ostride + j,
		     "output has been overwritten");

     fftw(validated_plan, 1, in2, 1, n, out2, 1, n);

     if (dir == FFTW_REAL_TO_COMPLEX) {
	  int n2 = (n + 1) / 2;
	  out3[0] = c_re(out2[0]);
	  for (i = 1; i < n2; ++i) {
	       out3[i] = c_re(out2[i]);
	       out3[n - i] = c_im(out2[i]);
	  }
	  if (n2 * 2 == n)
	       out3[n2] = c_re(out2[n2]);
     } else {
	  for (i = 0; i < n; ++i)
	       out3[i] = c_re(out2[i]);
     }

     for (j = 0; j < howmany; ++j)
	  CHECK(compute_error(out1 + j * n * ostride, ostride, out3, 1, n)
		< TOLERANCE,
		"test_out_of_place: wrong answer");
     WHEN_VERBOSE(2, printf("OK\n"));

     fftw_free(in1);
     fftw_free(in2);
     fftw_free(out1);
     fftw_free(out2);
     fftw_free(out3);
}

void test_in_place(int n, int istride,
		   int howmany, fftw_direction dir,
		   fftw_plan validated_plan, int specific)
{
     fftw_complex *in2, *out2;
     fftw_real *in1, *out1, *out3;
     fftw_plan plan;
     int i, j;
     int ostride = istride;
     int flags = measure_flag | wisdom_flag | FFTW_IN_PLACE;

     if (coinflip())
	  flags |= FFTW_THREADSAFE;

     in1 = (fftw_real *) fftw_malloc(istride * n * sizeof(fftw_real) * howmany);
     in2 = (fftw_complex *) fftw_malloc(n * sizeof(fftw_complex));
     out1 = in1;
     out2 = (fftw_complex *) fftw_malloc(n * sizeof(fftw_complex));
     out3 = (fftw_real *) fftw_malloc(n * sizeof(fftw_real));

     if (!specific)
	  plan = rfftw_create_plan(n, dir, flags);
     else
	  plan = rfftw_create_plan_specific(n, dir, flags,
					    in1, istride, out1, ostride);
     CHECK(plan != NULL, "can't create plan");

     /* generate random inputs */
     fill_random(in1, n, istride);
     for (j = 1; j < howmany; ++j)
	  for (i = 0; i < n; ++i)
	       in1[(j * n + i) * istride] = in1[i * istride];

     /* copy random inputs to complex array for comparison with fftw: */
     if (dir == FFTW_REAL_TO_COMPLEX)
	  for (i = 0; i < n; ++i) {
	       c_re(in2[i]) = in1[i * istride];
	       c_im(in2[i]) = 0.0;
     } else {
	  int n2 = (n + 1) / 2;
	  c_re(in2[0]) = in1[0];
	  c_im(in2[0]) = 0.0;
	  for (i = 1; i < n2; ++i) {
	       c_re(in2[i]) = in1[i * istride];
	       c_im(in2[i]) = in1[(n - i) * istride];
	  }
	  if (n2 * 2 == n) {
	       c_re(in2[n2]) = in1[n2 * istride];
	       c_im(in2[n2]) = 0.0;
	       ++i;
	  }
	  for (; i < n; ++i) {
	       c_re(in2[i]) = c_re(in2[n - i]);
	       c_im(in2[i]) = -c_im(in2[n - i]);
	  }
     }

     /* 
      * fill in other positions of the array, to make sure that
      * rfftw doesn't overwrite them 
      */
     for (j = 1; j < istride; ++j)
	  for (i = 0; i < n * howmany; ++i)
	       in1[i * istride + j] = i * istride + j;

     WHEN_VERBOSE(2, rfftw_print_plan(plan));

     /* fft-ize */
     if (howmany != 1 || istride != 1 || coinflip())
	  rfftw_threads(nthreads, plan, howmany, in1, istride, n * istride, 0, 0, 0);
     else
	  rfftw_threads_one(nthreads, plan, in1, NULL);

     rfftw_destroy_plan(plan);

     /* check for overwriting */
     for (j = 1; j < ostride; ++j)
	  for (i = 0; i < n * howmany; ++i)
	       CHECK(out1[i * ostride + j] == i * ostride + j,
		     "output has been overwritten");

     fftw(validated_plan, 1, in2, 1, n, out2, 1, n);

     if (dir == FFTW_REAL_TO_COMPLEX) {
	  int n2 = (n + 1) / 2;
	  out3[0] = c_re(out2[0]);
	  for (i = 1; i < n2; ++i) {
	       out3[i] = c_re(out2[i]);
	       out3[n - i] = c_im(out2[i]);
	  }
	  if (n2 * 2 == n)
	       out3[n2] = c_re(out2[n2]);
     } else {
	  for (i = 0; i < n; ++i)
	       out3[i] = c_re(out2[i]);
     }

     for (j = 0; j < howmany; ++j)
	  CHECK(compute_error(out1 + j * n * ostride, ostride, out3, 1, n)
		< TOLERANCE,
		"test_in_place: wrong answer");
     WHEN_VERBOSE(2, printf("OK\n"));

     fftw_free(in1);
     fftw_free(in2);
     fftw_free(out2);
     fftw_free(out3);
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

     /* produce a *good* plan (validated by Ergun's test procedure) */
     validated_plan_forward =
	 fftw_create_plan(n, FFTW_FORWARD, measure_flag | wisdom_flag);
     validated_plan_backward =
	 fftw_create_plan(n, FFTW_BACKWARD, measure_flag | wisdom_flag);
     CHECK(validated_plan_forward != NULL, "can't create plan");
     CHECK(validated_plan_backward != NULL, "can't create plan");

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

void testnd_out_of_place(int rank, int *n, fftwnd_plan validated_plan)
{
     int istride, ostride;
     int N, dim, i, j, k;
     int nc, nhc, nr;
     fftw_real *in1, *out3;
     fftw_complex *in2, *out1, *out2;
     fftwnd_plan p, ip;
     int flags = measure_flag | wisdom_flag;

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
     in1 = (fftw_real *) fftw_malloc(N * MAX_STRIDE * sizeof(fftw_real));
     out3 = (fftw_real *) fftw_malloc(N * MAX_STRIDE * sizeof(fftw_real));
     out1 = (fftw_complex *) fftw_malloc(nhc * nc * MAX_STRIDE
					 * sizeof(fftw_complex));
     in2 = (fftw_complex *) fftw_malloc(N * sizeof(fftw_complex));
     out2 = (fftw_complex *) fftw_malloc(N * sizeof(fftw_complex));

     p = rfftwnd_create_plan(rank, n, FFTW_REAL_TO_COMPLEX, flags);
     ip = rfftwnd_create_plan(rank, n, FFTW_COMPLEX_TO_REAL, flags);
     CHECK(p != NULL && ip != NULL, "can't create plan");

     for (istride = 1; istride <= MAX_STRIDE; ++istride) {
	  /* generate random inputs */
	  for (i = 0; i < nc; ++i)
	       for (j = 0; j < nr; ++j) {
		    c_re(in2[i * nr + j]) = DRAND();
		    c_im(in2[i * nr + j]) = 0.0;
		    for (k = 0; k < istride; ++k)
			 in1[(i * nr + j) * istride + k]
			     = c_re(in2[i * nr + j]);
	       }
	  for (i = 0; i < N * istride; ++i)
	       out3[i] = 0.0;

	  fftwnd(validated_plan, 1, in2, 1, 1, out2, 1, 1);

	  for (ostride = 1; ostride <= MAX_STRIDE; ++ostride) {
	       int howmany = (istride < ostride) ? istride : ostride;

	       WHEN_VERBOSE(2, printf("\n    testing stride %d/%d...",
				      istride, ostride));

	       if (howmany != 1 || istride != 1 || ostride != 1 || coinflip())
		    rfftwnd_threads_real_to_complex(nthreads, p, howmany, in1, istride, 1,
					    out1, ostride, 1);
	       else
		    rfftwnd_threads_one_real_to_complex(nthreads, p, in1, out1);

	       for (i = 0; i < nc; ++i)
		    for (k = 0; k < howmany; ++k)
			 CHECK(compute_error_complex(out1 + i * nhc * ostride + k,
						     ostride,
						     out2 + i * nr, 1,
						     nhc) < TOLERANCE,
			       "out-of-place (r2c): wrong answer");

	       if (howmany != 1 || istride != 1 || ostride != 1 || coinflip())
		    rfftwnd_threads_complex_to_real(nthreads, ip, howmany, out1, ostride, 1,
					    out3, istride, 1);
	       else
		    rfftwnd_threads_one_complex_to_real(nthreads, ip, out1, out3);

	       for (i = 0; i < N * istride; ++i)
		    out3[i] *= 1.0 / N;

	       if (istride == howmany)
		    CHECK(compute_error(out3, 1, in1, 1, N * istride)
			< TOLERANCE, "out-of-place (c2r): wrong answer");
	       for (i = 0; i < nc; ++i)
		    for (k = 0; k < howmany; ++k)
			 CHECK(compute_error(out3 + i * nr * istride + k,
					     istride,
					 (fftw_real *) (in2 + i * nr), 2,
					     nr) < TOLERANCE,
			   "out-of-place (c2r): wrong answer (check 2)");
	  }
     }

     rfftwnd_destroy_plan(p);
     rfftwnd_destroy_plan(ip);

     fftw_free(out3);
     fftw_free(out2);
     fftw_free(in2);
     fftw_free(out1);
     fftw_free(in1);
}

void testnd_in_place(int rank, int *n, fftwnd_plan validated_plan,
		     int alternate_api, int specific)
{
     int istride, ostride, howmany;
     int N, dim, i, j, k;
     int nc, nhc, nr;
     fftw_real *in1, *out3;
     fftw_complex *in2, *out1, *out2;
     fftwnd_plan p, ip;
     int flags = measure_flag | wisdom_flag | FFTW_IN_PLACE;

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
     in1 = (fftw_real *) fftw_malloc(2 * nhc * nc * MAX_STRIDE * sizeof(fftw_real));
     out3 = in1;
     out1 = (fftw_complex *) in1;
     in2 = (fftw_complex *) fftw_malloc(N * sizeof(fftw_complex));
     out2 = (fftw_complex *) fftw_malloc(N * sizeof(fftw_complex));

     if (alternate_api && specific && (rank == 2 || rank == 3)) {
	  if (rank == 2) {
	       p = rfftw2d_create_plan_specific(n[0], n[1],
					     FFTW_REAL_TO_COMPLEX, flags,
						in1, MAX_STRIDE, 0, 0);
	       ip = rfftw2d_create_plan_specific(n[0], n[1],
					     FFTW_COMPLEX_TO_REAL, flags,
						 in1, MAX_STRIDE, 0, 0);
	  } else {
	       p = rfftw3d_create_plan_specific(n[0], n[1], n[2],
					     FFTW_REAL_TO_COMPLEX, flags,
						in1, MAX_STRIDE, 0, 0);
	       ip = rfftw3d_create_plan_specific(n[0], n[1], n[2],
					     FFTW_COMPLEX_TO_REAL, flags,
						 in1, MAX_STRIDE, 0, 0);
	  }
     } else if (specific) {
	  p = rfftwnd_create_plan_specific(rank, n, FFTW_REAL_TO_COMPLEX,
					   flags,
				       in1, MAX_STRIDE, in1, MAX_STRIDE);
	  ip = rfftwnd_create_plan_specific(rank, n, FFTW_COMPLEX_TO_REAL,
					    flags,
				       in1, MAX_STRIDE, in1, MAX_STRIDE);
     } else if (alternate_api && (rank == 2 || rank == 3)) {
	  if (rank == 2) {
	       p = rfftw2d_create_plan(n[0], n[1], FFTW_REAL_TO_COMPLEX,
				       flags);
	       ip = rfftw2d_create_plan(n[0], n[1], FFTW_COMPLEX_TO_REAL,
					flags);
	  } else {
	       p = rfftw3d_create_plan(n[0], n[1], n[2], FFTW_REAL_TO_COMPLEX,
				       flags);
	       ip = rfftw3d_create_plan(n[0], n[1], n[2], FFTW_COMPLEX_TO_REAL,
					flags);
	  }
     } else {
	  p = rfftwnd_create_plan(rank, n, FFTW_REAL_TO_COMPLEX, flags);
	  ip = rfftwnd_create_plan(rank, n, FFTW_COMPLEX_TO_REAL, flags);
     }

     CHECK(p != NULL && ip != NULL, "can't create plan");

     for (i = 0; i < nc * nhc * 2 * MAX_STRIDE; ++i)
	  out3[i] = 0;

     for (istride = 1; istride <= MAX_STRIDE; ++istride) {
	  /* generate random inputs */
	  for (i = 0; i < nc; ++i)
	       for (j = 0; j < nr; ++j) {
		    c_re(in2[i * nr + j]) = DRAND();
		    c_im(in2[i * nr + j]) = 0.0;
		    for (k = 0; k < istride; ++k)
			 in1[(i * nhc * 2 + j) * istride + k]
			     = c_re(in2[i * nr + j]);
	       }

	  fftwnd(validated_plan, 1, in2, 1, 1, out2, 1, 1);

	  howmany = ostride = istride;

	  WHEN_VERBOSE(2, printf("\n    testing in-place stride %d...",
				 istride); fflush(stdout););

	  if (howmany != 1 || istride != 1 || ostride != 1 || coinflip())
	       rfftwnd_threads_real_to_complex(nthreads, p, howmany, in1, istride, 1,
				       out1, ostride, 1);
	  else
	       rfftwnd_threads_one_real_to_complex(nthreads, p, in1, NULL);

	  for (i = 0; i < nc; ++i)
	       for (k = 0; k < howmany; ++k)
		    CHECK(compute_error_complex(out1 + i * nhc * ostride + k,
						ostride,
						out2 + i * nr, 1,
						nhc) < TOLERANCE,
			  "in-place (r2c): wrong answer");

	  if (howmany != 1 || istride != 1 || ostride != 1 || coinflip())
	       rfftwnd_threads_complex_to_real(nthreads, ip, howmany, out1, ostride, 1,
				       out3, istride, 1);
	  else
	       rfftwnd_threads_one_complex_to_real(nthreads, ip, out1, NULL);

	  for (i = 0; i < nc * nhc * 2 * istride; ++i)
	       out3[i] *= 1.0 / N;

	  for (i = 0; i < nc; ++i)
	       for (k = 0; k < howmany; ++k)
		    CHECK(compute_error(out3 + i * nhc * 2 * istride + k,
					istride,
					(fftw_real *) (in2 + i * nr), 2,
					nr) < TOLERANCE,
			  "in-place (c2r): wrong answer (check 2)");
     }

     rfftwnd_destroy_plan(p);
     rfftwnd_destroy_plan(ip);

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

     testnd_out_of_place(sz.rank, sz.narray, validated_plan);
     testnd_in_place(sz.rank, sz.narray,
		     validated_plan, alt_api, specific);

     fftwnd_destroy_plan(validated_plan);
}

/*************************************************
 * planner tests
 *************************************************/

void test_planner(int rank)
{
     WHEN_VERBOSE(1, printf("Use rfftw_test to test the planner.\n"););
}

/*************************************************
 * test initialization
 *************************************************/

void test_init(int *argc, char **argv)
{
     int i;

     /*     if (*argc >= 2)
          nthreads = atoi(argv[1]);

     if (nthreads <= 0) {
          fprintf(stderr, "Usage: fftw_threads_test <nthreads> [ options ]\n");
          exit(EXIT_FAILURE);
     }
     for (i = 2; i < *argc; ++i)
          argv[i - 1] = argv[i];
     *argc -= 1;
     */
     nthreads=4;
printf("using 4 threads!\n");
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
