/*
 * Copyright (c) 1997 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to use, copy, modify, and distribute the Software without
 * restriction, provided the Software, including any modified copies made
 * under this license, is not distributed for a fee, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE MASSACHUSETTS INSTITUTE OF TECHNOLOGY BE LIABLE
 * FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 * Except as contained in this notice, the name of the Massachusetts
 * Institute of Technology shall not be used in advertising or otherwise
 * to promote the sale, use or other dealings in this Software without
 * prior written authorization from the Massachusetts Institute of
 * Technology.
 *  
 */

/*
 * timer.c -- this file measures the execution time of 
 *            ffts.  This information is used by the planner.
 */

/* $Id$ */

#include <time.h>
#include <fftw.h>
#include <math.h>
#include <stdlib.h>

/*
 * The timer keeps doubling the number of iterations
 * until the program runs for more than FFTW_TIME_MIN
 */
double fftw_measure_runtime(fftw_plan plan)
{
     FFTW_COMPLEX *in, *out;
     fftw_time begin, end;
     double t;
     int i, iter;
     int n;

     n = plan->n;

     iter = 1;

retry:
     in = (FFTW_COMPLEX *) fftw_malloc(n * sizeof(FFTW_COMPLEX));
     out = (FFTW_COMPLEX *) fftw_malloc(n * sizeof(FFTW_COMPLEX));

     begin = fftw_get_time();
     for (i = 0; i < iter; ++i) {
	  int j;

	  /* generate random inputs */
	  for (j = 0; j < n; ++j) {
	       c_re(in[j]) = 1.0;
	       c_im(in[j]) = 32.432;
	  }

	  fftw(plan, 1, in, 1, 0, out, 1, 0);
     }
     end = fftw_get_time();

     t = fftw_time_to_sec(fftw_time_diff(end,begin));

     fftw_free(in);
     fftw_free(out);

     if (t < FFTW_TIME_MIN) {
	  iter *= 2;
	  /* 
	   * See D. E. Knuth, Structured Programming with GOTO Statements,
	   * Computing Surveys (6), December 1974, for a justification
	   * of this `goto' in the `n + 1/2' loop.
	   */
	  goto retry;
     }

     return t / (double)iter;
}

#if defined(MAC) || defined(macintosh)

/* Use Macintosh Time Manager to get the time: */

#pragma only_std_keywords off  /* make sure compiler (CW) recognizes the pascal
				  keywords that are in Timer.h */

#include <Timer.h>

#pragma only_std_keywords reset

fftw_time get_Mac_microseconds(void)
{
     fftw_time t;
     UnsignedWide microsec;	/* 
				 * microsec.lo and microsec.hi are
				 * unsigned long's, and are the two parts
				 * of a 64 bit unsigned integer 
				 */

     Microseconds(&microsec);	/* get time in microseconds */

     /* store lo and hi words into our structure: */
     t.lo = microsec.lo; t.hi = microsec.hi;

     return t;
}

fftw_time fftw_time_diff(fftw_time t1, fftw_time t2)
/* This function takes the difference t1 - t2 of two 64 bit
   integers, represented by the 32 bit lo and hi words.
   if t1 < t2, returns 0. */
{
     fftw_time diff;

     if (t1.hi < t2.hi) { /* something is wrong...t1 < t2! */
	  diff.hi = diff.lo = 0;
	  return diff;
     }
     else
	  diff.hi = t1.hi - t2.hi;

     if (t1.lo < t2.lo) {
	  if (diff.hi > 0)
	       diff.hi -= 1; /* carry */
	  else { /* something is wrong...t1 < t2! */
	       diff.hi = diff.lo = 0;
	       return diff;
	  }
     }
     
     diff.lo = t1.lo - t2.lo;

     return diff;
}

#endif

#ifdef __WIN32__
#include <windows.h>

static LARGE_INTEGER gFreq;
static int gHaveHiResTimer = 0;
static int gFirstTime = 1;

unsigned long GetPerfTime(void)
{
     LARGE_INTEGER lCounter;

     if (gFirstTime) {
	  gFirstTime = 0;

	  if (QueryPerformanceFrequency(&gFreq)) {
	       gHaveHiResTimer = 1;
	  }
     }
     if (gHaveHiResTimer) {
	  QueryPerformanceCounter(&lCounter);
	  return lCounter.u.LowPart;
     } else {
	  return (unsigned long) clock();
     }
}

double GetPerfSec(double pTime)
{
     if (gHaveHiResTimer) {
	  return pTime / gFreq.u.LowPart;	// assumes HighPart==0

     } else {
	  return pTime / CLOCKS_PER_SEC;
     }
}

#endif				/* __WIN32__ */
