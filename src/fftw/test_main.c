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
 * test_main.c: driver for test programs (linked with fftw_test.c/rfftw_test.c)
 */

/* $Id$ */
#include <fftw-int.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "test_main.h"

#ifdef HAVE_GETOPT
#  ifdef HAVE_GETOPT_H
#    include <getopt.h>
#  elif defined(HAVE_UNISTD_H)
#    include <unistd.h>
#  endif
#endif
#include <unistd.h>
#define HAVE_GETOPT
double mydrand(void)
{
     double d = rand();
     return (d / (double) RAND_MAX) - .5;
}

/* return random 0 or non-zero */
int coinflip(void)
{
     return (rand() & 8192);	/* higher-order bits are often more random 
				 * 
				 */
}

/* parse a string of the form N1xN2x... and return a size structure */
struct size parse_size(char *s)
{
     struct size sz;
     int n;

     sz.rank = 0;
     sz.is_nd = 0;

     if (*s == 'x' || *s == 'X' || *s == '*') {
	  ++s;
	  /* force ND transform of rank 1 */
	  sz.is_nd = 1;
     }
   accept_digit:
     n = 0;

     CHECK(isdigit(*s),
	   "invalid digit");

     while (isdigit(*s)) {
	  n = n * 10 + (*s - '0');
	  ++s;
     }

     sz.narray[sz.rank] = n;
     ++sz.rank;

     CHECK(sz.rank < MAX_CMDLINE_RANK,
	   "maximum rank exceeded");

     if (*s == 'x' || *s == 'X' || *s == '*') {
	  ++s;
	  goto accept_digit;
     }
     /* force ND transform if rank > 1 */
     if (sz.rank > 1)
	  sz.is_nd = 1;
     return sz;
}

/*******************
 * global variables
 *******************/
int verbose;
int wisdom_flag, measure_flag;
int speed_flag = FFTW_MEASURE;
int chk_mem_leak;
int paranoid;
int howmany_fields = 1;
int max_iterations = 0;		/* 
				 * maximum number of iterations to perform
				 * in "infinite" tests--default (0) means
				 * no limit 
				 */

/* When testing MPI stuff, only one process gets to do I/O: */
int io_okay = 1;
#define my_printf if (io_okay) printf
#define my_fprintf if (io_okay) fprintf
#define my_fflush if (io_okay) fflush

/*******************
 * procedures
 *******************/

/* smart time printer */
char *smart_sprint_time(double x)
{
     static char buf[128];

     if (x < 1.0E-6)
	  sprintf(buf, "%f ns", x * 1.0E9);
     else if (x < 1.0E-3)
	  sprintf(buf, "%f us", x * 1.0E6);
     else if (x < 1.0)
	  sprintf(buf, "%f ms", x * 1.0E3);
     else
	  sprintf(buf, "%f s", x);

     return buf;
}

/* greet the user */
/* jokes stolen from http://whereis.mit.edu/bin/map */
void please_wait(void)
{
     int i;
     const char *s[] =
     {
       "(while a large software vendor in Seattle takes over the world)",
       "(and remember, this is faster than Java)",
       "(and dream of faster computers)",
       "(checking the gravitational constant in your locale)",
       "(at least you are not on hold)",
       "(while X11 grows by another kilobyte)",
       "(while Windows NT reboots)",
       "(correcting for the phase of the moon)",
       "(your call is important to us)",
       "(while the Linux user-base doubles)",
       "(while you decide where you want to go tomorrow)",
       "(exorcising evil spirits)",
     };
     int choices = sizeof(s) / sizeof(*s);

     i = rand() % choices;
     my_printf("Please wait %s.\n", s[i < 0 ? -i : i]);
}

void please_wait_forever(void)
{
     int i;
     const char *s[] =
     {
	  "(but it won't crash, either)",
	  "(at least in theory)",
	  "(please be patient)",
	  "(our next release will complete it more quickly)",
#if defined(__WIN32__) || defined(WIN32) || defined(_WINDOWS)
	  "(by the way, Linux executes infinite loops faster)",
#endif
     };
     int choices = sizeof(s) / sizeof(*s);

     if (!max_iterations) {
	  i = rand() % choices;
	  my_printf("This test does not terminate %s.\n", s[i < 0 ? -i : i]);
     } else {
	  my_printf("This test will run for %d iterations.\n", max_iterations);
	  please_wait();
     }
}

/*************************************************
 * Speed tests
 *************************************************/

double mflops(double t, int N)
{
     return (5.0 * N * log((double) N) / (log(2.0) * t * 1.0e6));
}

void print_dims(struct size sz)
{
     int i;

     my_printf("%d", sz.narray[0]);
     for (i = 1; i < sz.rank; ++i)
	  my_printf("x%d", sz.narray[i]);
}

void test_speed(int n)
{
     int specific;

     please_wait();

     if (howmany_fields > 1)
	  WHEN_VERBOSE(1, my_printf("TIMING MULTIPLE-FIELD FFT: "
				 "howmany=%d, stride=%d, dist=%d\n\n",
				 howmany_fields, howmany_fields, 1));

     for (specific = 0; specific <= 1; ++specific) {
	  WHEN_VERBOSE(1,
	   my_printf("SPEED TEST: n = %d, FFTW_FORWARD, out of place, %s\n",
		  n, SPECIFICP(specific)));
	  test_speed_aux(n, FFTW_FORWARD, 0, specific);

	  WHEN_VERBOSE(1,
	       my_printf("SPEED TEST: n = %d, FFTW_FORWARD, in place, %s\n",
		      n, SPECIFICP(specific)));
	  test_speed_aux(n, FFTW_FORWARD, FFTW_IN_PLACE, specific);

	  WHEN_VERBOSE(1,
	  my_printf("SPEED TEST: n = %d, FFTW_BACKWARD, out of place, %s\n",
		 n, SPECIFICP(specific)));
	  test_speed_aux(n, FFTW_BACKWARD, 0, specific);

	  WHEN_VERBOSE(1,
	      my_printf("SPEED TEST: n = %d, FFTW_BACKWARD, in place, %s\n",
		     n, SPECIFICP(specific)));
	  test_speed_aux(n, FFTW_BACKWARD, FFTW_IN_PLACE, specific);
     }
}

void test_speed_nd(struct size sz)
{
     int specific;

     please_wait();

     if (howmany_fields > 1)
	  WHEN_VERBOSE(1, my_printf("TIMING MULTIPLE-FIELD FFT: "
				 "howmany=%d, stride=%d, dist=%d\n\n",
				 howmany_fields, howmany_fields, 1));

     for (specific = 0; specific <= 1; ++specific) {
	  my_printf("SPEED TEST: ");
	  WHEN_VERBOSE(1, print_dims(sz));
	  WHEN_VERBOSE(1, my_printf(", FFTW_FORWARD, in place, %s\n",
				 SPECIFICP(specific)));
	  test_speed_nd_aux(sz, FFTW_FORWARD,
			    FFTW_IN_PLACE, specific);

	  WHEN_VERBOSE(1, my_printf("SPEED TEST: "));
	  print_dims(sz);
	  WHEN_VERBOSE(1, my_printf(", FFTW_BACKWARD, in place, %s\n",
				 SPECIFICP(specific)));
	  test_speed_nd_aux(sz, FFTW_BACKWARD, FFTW_IN_PLACE, specific);
     }
}

/*************************************************
 * correctness tests
 *************************************************/

double compute_error_complex(fftw_complex * A, int astride,
			     fftw_complex * B, int bstride, int n)
{
     /* compute the relative error */
     double error = 0.0;
     int i;

     for (i = 0; i < n; ++i) {
	  double a;
	  double mag;
	  a = sqrt(SQR(c_re(A[i * astride]) - c_re(B[i * bstride])) +
		   SQR(c_im(A[i * astride]) - c_im(B[i * bstride])));
	  mag = 0.5 * (sqrt(SQR(c_re(A[i * astride]))
			    + SQR(c_im(A[i * astride]))) +
		       sqrt(SQR(c_re(B[i * bstride]))
			    + SQR(c_im(B[i * bstride])))) + TOLERANCE;

	  a /= mag;
	  if (a > error)
	       error = a;

#ifdef HAVE_ISNAN
	  CHECK(!isnan(a), "NaN in answer");
#endif
     }
     return error;
}

/* test forever */
void test_all(void)
{
     int n;

     please_wait_forever();
     for (n = 1; !max_iterations || n <= max_iterations; ++n) {
	  test_correctness(n);
	  if (!(wisdom_flag & FFTW_USE_WISDOM) && chk_mem_leak)
	       fftw_check_memory_leaks();
     }
}

#define MAX_FACTOR 13

int rand_small_factors(int N)
{
     int f, n = 1;

     f = rand() % MAX_FACTOR + 1;

     while (n * f <= N) {
	  n *= f;
	  f = rand() % MAX_FACTOR + 1;
     }

     return n;
}

#define MAX_N 16384

struct size random_dims(int rank)
{
     int maxsize, dim;
     double maxsize_d;
     struct size sz;

     /* workaround to weird gcc warning */
     maxsize_d = pow((double) (rank == 1 ? MAX_N / 4 : MAX_N),
		     1.0 / (double) rank);
     maxsize = (int) maxsize_d;

     if (maxsize < 1)
	  maxsize = 1;

     sz.rank = rank;
     for (dim = 0; dim < rank; ++dim)
	  sz.narray[dim] = rand_small_factors(maxsize);

     return sz;
}

void test_random(void)
{
     static int counter = 0;
     struct size sz;

     if ((++counter) % 16 == 0) {
	  sz.rank = 1;
	  sz.narray[0] = rand() % (MAX_N / 16) + 1;
     } else {
	  sz = random_dims(1);
     }

     test_correctness(sz.narray[0]);
}

/*************************************************
 * multi-dimensional correctness tests
 *************************************************/

void testnd_correctness_both(struct size sz,
			     int alt_api, int specific, int force_buf)
{
     WHEN_VERBOSE(1,
		  my_printf("Testing nd correctness for size = ");
		  print_dims(sz);
		  my_printf("...");
		  my_fflush(stdout));

     if (alt_api)
	  WHEN_VERBOSE(1, my_printf("alt. api..."));
     if (specific)
	  WHEN_VERBOSE(1, my_printf("specific..."));
     if (force_buf)
	  WHEN_VERBOSE(1, my_printf("force buf..."));

     testnd_correctness(sz, FFTW_FORWARD, alt_api, specific, force_buf);
     testnd_correctness(sz, FFTW_BACKWARD, alt_api, specific, force_buf);

     WHEN_VERBOSE(1, my_printf("OK\n"));
}

void testnd_correctness_aux(struct size sz)
{
     int alt_api, specific, force_buf;

     for (alt_api = 0; alt_api <= 1; ++alt_api)
	  for (specific = 0; specific <= 1; ++specific)
	       for (force_buf = 0; force_buf <= 1; ++force_buf)
		    testnd_correctness_both(sz, alt_api, specific,
					    force_buf);
}

void testnd_correctness_square(int rank, int size)
{
     struct size sz;
     int alt_api, specific, force_buf;
     int i;

     sz.rank = rank;
     for (i = 0; i < rank; ++i)
	  sz.narray[i] = size;

     for (alt_api = 0; alt_api <= 1; ++alt_api)
	  for (specific = 0; specific <= 1; ++specific)
	       for (force_buf = 0; force_buf <= 1; ++force_buf)
		    testnd_correctness_both(sz, alt_api, 
					    specific, force_buf);

}

void testnd_random(int rank)
{
     struct size sz;

     sz = random_dims(rank);
     testnd_correctness_both(sz, coinflip(), coinflip(), coinflip());
}

/* loop forever */
void test_all_random(int rank)
{
     int counter;
     please_wait_forever();

     for (counter = 0; !max_iterations || counter < max_iterations; ++counter) {
	  if (rank > 0)
	       testnd_random(rank);
	  else if ((counter) % 2 == 0)
	       test_random();
	  else
	       testnd_random(rand() % MAX_RANK + 1);
     }

}

int pow2sqrt(int n)
/* return greatest power of two <= sqrt(n) */
{
     int s = 1;

     while (s * s * 4 <= n)
	  s *= 2;
     return s;
}

/* test forever */
void testnd_all(int rank)
{
     int n;

     please_wait_forever();

     for (n = 1; !max_iterations || n <= max_iterations; ++n)
	  testnd_correctness_square(rank, n);
}

fftw_direction random_dir(void)
{
     if (coinflip())
	  return FFTW_FORWARD;
     else
	  return FFTW_BACKWARD;
}

/*************************************************
 * timer tests
 *************************************************/

static int hack_sum_i;

void negative_time(void)
{
     my_fprintf(stderr,
	     "* PROBLEM: I measured a negative time interval.\n"
	     "* Please make sure you defined the timer correctly\n"
	     "* or contact fftw@theory.lcs.mit.edu for help.\n");
}

/*
 * paranoid test to see if time is monotonic.  If not, you are
 * really in trouble
 */
void test_timer_paranoid(void)
{
     fftw_time start_t, end_t;
     double sec;
     int i;

     start_t = fftw_get_time();

     /* waste some time */
     for (i = 0; i < 10000; ++i)
	  hack_sum_i = i;

     end_t = fftw_get_time();
     sec = fftw_time_to_sec(fftw_time_diff(end_t, start_t));
     if (sec < 0.0)
	  negative_time();
}

/* compute useful numbers */
static int fib(int n)
{
     if (n < 2)
	  return n;
     else {
	  int x, y;
	  x = fib(n - 1);
	  y = fib(n - 2);
	  return x + y;
     }
}

static int hack_fib;

void test_timer(void)
{
     double times[32], acc, min_time = 10000.00;
     unsigned long iters, iter;
     fftw_time begin, end, start;
     double t, tmax, tmin;
     int last = 0, i, repeat;

     please_wait();
     test_timer_paranoid();

     start = fftw_get_time();

     for (i = 0; i < 32; i++) {
	  iters = 1 << i;
	  tmin = 1.0E10;
	  tmax = -1.0E10;

	  for (repeat = 0; repeat < FFTW_TIME_REPEAT; ++repeat) {
	       begin = fftw_get_time();
	       for (iter = 0; iter < iters; ++iter) {
		    hack_fib = fib(10);
	       }
	       end = fftw_get_time();

	       t = fftw_time_to_sec(fftw_time_diff(end, begin));
	       if (t < tmin)
		    tmin = t;
	       if (t > tmax)
		    tmax = t;

	       /* do not run for too long */
	       t = fftw_time_to_sec(fftw_time_diff(end, start));
	       if (t > FFTW_TIME_LIMIT)
		    break;
	  }

	  if (tmin < 0.0)
	       negative_time();

	  times[i] = tmin;

	  WHEN_VERBOSE(2,
		  my_printf("Number of iterations = 2^%d = %lu, time = %g, "
			 "time/iter = %g\n",
			 i, iters, times[i],
			 times[i] / iters));
	  WHEN_VERBOSE(2,
		   my_printf("   (out of %d tries, tmin = %g, tmax = %g)\n",
			  FFTW_TIME_REPEAT, tmin, tmax));

	  last = i;
	  if (times[i] > 10.0)
	       break;
     }

     /* 
      * at this point, `last' is the last valid element in the
      * `times' array.
      */

     for (i = 0; i <= last; ++i)
	  if (times[i] > 0.0 && times[i] < min_time)
	       min_time = times[i];

     WHEN_VERBOSE(1, my_printf("\nMinimum resolvable time interval = %g seconds.\n\n",
			    min_time));

     for (acc = 0.1; acc > 0.0005; acc *= 0.1) {
	  double t_final;
	  t_final = times[last] / (1 << last);

	  for (i = last; i >= 0; --i) {
	       double t_cur, error;
	       iters = 1 << i;
	       t_cur = times[i] / iters;
	       error = (t_cur - t_final) / t_final;
	       if (error < 0.0)
		    error = -error;
	       if (error > acc)
		    break;
	  }

	  ++i;

	  WHEN_VERBOSE(1,
	      my_printf("Minimum time for %g%% consistency = %g seconds.\n",
		     acc * 100.0, times[i]));
     }
     WHEN_VERBOSE(1,
	      my_printf("\nMinimum time used in FFTW timing (FFTW_TIME_MIN)"
		     " = %g seconds.\n", FFTW_TIME_MIN));
}

/*************************************************
 * help
 *************************************************/

#ifdef HAVE_GETOPT_LONG
#  define WHEN_LONG_OPTIONS(x) x
#else
#  define WHEN_LONG_OPTIONS(x) ""
#endif

static void usage(int exit_when_done)
{
     my_printf("Usage:  %s_test [options]\n", fftw_prefix);
     my_printf(WHEN_LONG_OPTIONS("  --speed=<n>\n")
	       "  -s <n>    : test speed for size n\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --correctness=<n>\n")
	       "  -c <n>    : test correctness for size n\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --random=<rank>>\n")
	       "  -r <rank> : test correctness for random sizes "
	       "(does not terminate)\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --all=<rank>\n")
	       "  -a <rank> : test correctness for all sizes "
	       "(does not terminate)\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --fields=<n>\n")
	       "  -f <n>    : n fields ('howmany' param) in speed tests\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --planner=<rank>\n")
	       "  -p <rank> : test planner\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --measure\n")
	       "  -m        : use FFTW_MEASURE in correctness tests\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --estimate\n")
	       "  -e        : use FFTW_ESTIMATE in speed tests\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --wisdom=<file>\n")
	       "  -w <file> : use wisdom & read/write it from/to file\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --timer\n")
	       "  -t        : test timer resolution\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --x-repeat=<n>\n")
	       "  -x <n>    : run non-terminating tests (-r, -a) only n times\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --paranoid\n")
	       "  -P        : enable paranoid tests\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --verbose\n")
	       "  -v        : verbose output for subsequent options\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --version\n")
	       "  -V        : print FFTW version information\n");
     my_printf(WHEN_LONG_OPTIONS("\n  --help\n")
	       "  -h        : this help\n");
#ifndef HAVE_GETOPT
     my_printf("(When run with no arguments, an interactive mode is used.)\n");
#endif
     if (exit_when_done)
	  exit(EXIT_FAILURE);
}

char wfname[128];

void handle_option(char opt, char *optarg)
{
     FILE *wf;
     struct size sz;
     int rank, n;

     switch (opt) {
	 case 's':
	      sz = parse_size(optarg);
	      if (!sz.is_nd)
		   test_speed(sz.narray[0]);
	      else
		   test_speed_nd(sz);
	      break;

	 case 'c':
	      sz = parse_size(optarg);
	      if (!sz.is_nd)
		   test_correctness(sz.narray[0]);
	      else
		   testnd_correctness_aux(sz);
	      break;

	 case 'p':
	      rank = atoi(optarg);
	      test_planner(rank);
	      break;

	 case 'P':
	      paranoid = 1;
	      enter_paranoid_mode();
	      break;

	 case 'r':
	      rank = atoi(optarg);
	      test_all_random(rank);
	      break;

	 case 'a':
	      rank = atoi(optarg);
	      if (rank == 0)
		   test_all();
	      else
		   testnd_all(rank);
	      break;

	 case 't':
	      test_timer();
	      break;

	 case 'f':
	      n = atoi(optarg);
	      CHECK(n > 0, "-f requires a positive integer argument");
	      howmany_fields = n;
	      break;

	 case 'm':
	      measure_flag = FFTW_MEASURE;
	      break;

	 case 'e':
	      speed_flag = FFTW_ESTIMATE;
	      break;

	 case 'w':
	      wisdom_flag = FFTW_USE_WISDOM;
	      strcpy(wfname, optarg);
	      wf = fopen(wfname, "r");
	      if (wf == 0) {
		   my_printf("Couldn't open wisdom file \"%s\".\n", wfname);
		   my_printf("This file will be created upon completion.\n");
	      } else {
		   CHECK(FFTW_SUCCESS == fftw_import_wisdom_from_file(wf),
			 "invalid wisdom file format");
		   fclose(wf);
	      }
	      break;

	 case 'v':
	      verbose++;
	      break;

	 case 'V':
	      my_printf("%s\n", fftw_version);
	      my_printf("%s test program, compiled in %s precision.\n",
			fftw_prefix,
			sizeof(fftw_real) == sizeof(double) ? "double"
			: (sizeof(fftw_real) == sizeof(float) ? "single"
			   : "unknown"));
	      my_printf(
  "\nCopyright (C) Massachusetts Institute of Technology.\n"
  "FFTW comes with ABSOLUTELY NO WARRANTY.  This is free software, and\n"
  "you are welcome to redistribute it under the terms of the GNU\n"
  "General Public License.  For more information, see the file COPYING or\n"
  "the GNU web site at http://www.gnu.org.\n"
  "\nFor more information regarding FFTW, or to download the latest version,\n"
  "see the FFTW home page at http://theory.lcs.mit.edu/~fftw.\n");

	      break;

	 case 'x':
	      n = atoi(optarg);
	      CHECK(n > 0, "-x requires a positive integer argument");
	      max_iterations = n;
	      break;

	 case 'h':
	      usage(FALSE);
	      break;

	 default:
	      usage(TRUE);
     }

     /* every test must free all the used FFTW memory */
     if (!(wisdom_flag & FFTW_USE_WISDOM) && chk_mem_leak)
	  fftw_check_memory_leaks();
}



short askuser(const char *s)
{
     char line[200] = "", c;
     int i, count = 0;

     do {
	  if (count++ > 0)
	       my_printf("Invalid response.  Please enter \"y\" or \"n\".\n");
	  my_printf("%s (y/n) ", s);
	  /* skip blank lines */
	  while (line[0] == 0 || line[0] == '\n')
	       fgets(line, 200, stdin);
	  for (i = 0; line[i] && (line[i] == ' ' || line[i] == '\t'); ++i);
	  c = line[i];
     } while (c != 'n' && c != 'N' && c != 'y' && c != 'Y');

     return (c == 'y' || c == 'Y');
}

/* Standard function to get the next command-line argument for the program.
   Returns the option character (or -1 if there are no more options),
   and the option argument (if any) in argval, which is an array of length
   at least argval_maxlen.

   The test programs need to implement a function get_option with the
   same arguments as this one, which will typically just call
   default_get_option.

   The reason we need to put this in a separate function is that the MPI
   test programs can't rely on all of the processes having direct access
   to the program arguments--they need to pass them as explicit messages
   from the master process.   Sigh. */
int default_get_option(int argc, char **argv, char *argval, int argval_maxlen)
{
     int c = -1;

     if (argc <= 1)
	  usage(TRUE);

#ifdef HAVE_GETOPT
     {
	  const char short_options[] = "s:c:w:f:p:Pa:r:tvVmehx:";
	  extern char *optarg;
	  extern int optind;
	  
#  if defined(HAVE_GETOPT_LONG) && defined(HAVE_GETOPT_H)
	  {
	       int option_index;
	       const struct option long_options[] = {
		    {"speed", 1, 0, 's'},
		    {"correctness", 1, 0, 'c'},
		    {"wisdom", 1, 0, 'w'},
		    {"fields", 1, 0, 'f'},
		    {"planner", 1, 0, 'p'},
		    {"paranoid", 0, 0, 'P'},
		    {"all", 1, 0, 'a'},
		    {"random", 1, 0, 'r'},
		    {"timer", 0, 0, 't'},
		    {"verbose", 0, 0, 'v'},
		    {"version", 0, 0, 'V'},
		    {"measure", 0, 0, 'm'},
		    {"estimate", 0, 0, 'e'},
		    {"help", 0, 0, 'h'},
		    {"x-repeat", 1, 0, 'x'},
		    {0, 0, 0, 0}
	       };
	       
	       c = getopt_long(argc, argv, short_options, long_options,
			       &option_index);
	  }
#  else /* not HAVE_GETOPT_LONG */
	  c = getopt(argc, argv, short_options);
#  endif /* not HAVE_GETOPT_LONG */
	  
	  if (c == -1 && argc != optind)
	       usage(TRUE);  /* there were invalid args; print usage info */

	  if (optarg) {
	       strncpy(argval, optarg, argval_maxlen - 1);
	       argval[argval_maxlen - 1] = 0;
	  }
	  else
	       argval[0] = 0;
     }
#endif /* HAVE_GETOPT */

     return c;
}

int main(int argc, char *argv[])
{
     verbose = 1;
     wisdom_flag = 0;
     measure_flag = FFTW_ESTIMATE;
     chk_mem_leak = 1;
     paranoid = 0;

#ifdef DETERMINISTIC
     srand(1123);
#else
     srand((unsigned int) time(NULL));
#endif

     test_init(&argc, argv);

     /* 
      * To parse the command line, we use getopt, but this does not seem
      * to be in the ANSI standard (it is only available on UNIX,
      * apparently). 
      */
#ifndef HAVE_GETOPT
     if (argc > 1)
	  my_printf("Sorry, command-line arguments are not available on\n"
		 "this system.  Run fftw_test with no arguments to\n"
		 "use it in interactive mode.\n");

     if (argc <= 1) {
	  int n = 0;
	  char s[128] = "";

	  usage(FALSE);

	  my_printf("\n");

	  if (askuser("Perform random correctness tests (non-terminating)?"))
	       handle_option('r', "0");

	  if (askuser("Verbose output?"))
	       handle_option('v', "");
	  if (askuser("Paranoid test?"))
	       handle_option('P', "");

	  if (askuser("Use/test wisdom?")) {
	       my_printf("  Enter wisdom file name to use: ");
	       fgets(s, 128, stdin);
	       handle_option('w', s);
	  }
	  if (askuser("Test correctness?")) {
	       if (askuser("  -- for all sizes?")) 
		    handle_option('a', "");
	       else {
		    my_printf("  Enter n: ");
		    fgets(s, 128, stdin);
		    handle_option('c', s);
	       }
	  }
	  if (askuser("Test speed?")) {
	       my_printf("  Enter n: ");
	       fgets(s, 128, stdin);
	       handle_option('s', s);
	  }
	  if (askuser("Test planner?"))
	       handle_option('p', "");
	  if (askuser("Test timer?"))
	       handle_option('t', "");
     }
#else				/* 
				 * read command-line args using getopt 
				 * facility  
				 */
     {
	  char option_arg[128];
	  int c;

	  while ((c = get_option(argc, argv, option_arg, 128)) != -1)
	       handle_option(c, option_arg);
     }
#endif

     if (wisdom_flag & FFTW_USE_WISDOM) {
	  char *ws;
	  FILE *wf;

	  ws = fftw_export_wisdom_to_string();
	  CHECK(ws != 0, "error exporting wisdom to string");
	  my_printf("\nAccumulated wisdom:\n     %s\n", ws);
	  fftw_forget_wisdom();
	  CHECK(FFTW_SUCCESS == fftw_import_wisdom_from_string(ws),
		"unexpected error reading in wisdom from string");
	  fftw_free(ws);

	  if (io_okay) {
	       wf = fopen(wfname, "w");
	       CHECK(wf != 0, "error creating wisdom file");
	       fftw_export_wisdom_to_file(wf);
	       fclose(wf);
	  }
     }
     /* make sure to dispose of wisdom before checking for memory leaks */
     fftw_forget_wisdom();

     fftw_check_memory_leaks();
     if (io_okay)
	  fftw_print_max_memory_usage();

     test_finish();

     return EXIT_SUCCESS;
}
