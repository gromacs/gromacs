
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

/*
 * timer.c -- this file measures the execution time of 
 *            ffts.  This information is used by the planner.
 */

/* $Id$ */

#include <time.h>
#include <fftw-int.h>
#include <math.h>
#include <stdlib.h>

/********************* System-specific Timing Support *********************/

#if defined(HAVE_MAC_TIMER) && !defined(HAVE_MAC_PCI_TIMER)

/* Use Macintosh Time Manager to get the time: */

/*
 * make sure compiler (CW) recognizes the pascal keywords that are in
 * Timer.h
 */
#pragma only_std_keywords off	

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
     t.lo = microsec.lo;
     t.hi = microsec.hi;

     return t;
}

fftw_time fftw_time_diff(fftw_time t1, fftw_time t2)
/*
 * This function takes the difference t1 - t2 of two 64 bit
 * integers, represented by the 32 bit lo and hi words.
 * if t1 < t2, returns 0. 
 */
{
     fftw_time diff;

     if (t1.hi < t2.hi) {	/* something is wrong...t1 < t2! */
	  diff.hi = diff.lo = 0;
	  return diff;
     } else
	  diff.hi = t1.hi - t2.hi;

     if (t1.lo < t2.lo) {
	  if (diff.hi > 0)
	       diff.hi -= 1;	/* carry */
	  else {		/* something is wrong...t1 < t2! */
	       diff.hi = diff.lo = 0;
	       return diff;
	  }
     }
     diff.lo = t1.lo - t2.lo;

     return diff;
}

#endif

#ifdef HAVE_WIN32_TIMER
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

#endif				/* HAVE_WIN32_TIMER */

#if defined(FFTW_USE_GETTIMEOFDAY)

/* timer support routines for systems having gettimeofday */

#ifdef HAVE_BSDGETTIMEOFDAY
#define gettimeofday BSDgettimeofday
#endif

fftw_time fftw_gettimeofday_get_time(void)
{
     struct timeval tv;
     gettimeofday(&tv, 0);
     return tv;
}

fftw_time fftw_gettimeofday_time_diff(fftw_time t1, fftw_time t2)
{
     fftw_time diff;

     diff.tv_sec = t1.tv_sec - t2.tv_sec;
     diff.tv_usec = t1.tv_usec - t2.tv_usec;
     /* normalize */
     while (diff.tv_usec < 0) {
	  diff.tv_usec += 1000000L;
	  diff.tv_sec -= 1;
     }

     return diff;
}
#endif
