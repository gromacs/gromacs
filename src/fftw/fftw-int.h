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

/* fftw.h -- system-wide definitions */
/* $Id$ */

#ifndef FFTW_INT_H
#define FFTW_INT_H
#include "fftw.h"

#ifdef __cplusplus
extern "C" {
#else
#endif				/* __cplusplus */

/****************************************************************************/
/*                            Private Functions                             */
/****************************************************************************/

extern fftw_twiddle *fftw_create_twiddle(int n, const fftw_codelet_desc *d);
extern void fftw_destroy_twiddle(fftw_twiddle *tw);

extern void fftw_strided_copy(int, fftw_complex *, int, fftw_complex *);
extern void fftw_executor_simple(int, const fftw_complex *, fftw_complex *,
				 fftw_plan_node *, int, int);

extern fftwnd_plan fftwnd_create_plan_aux(int rank, const int *n,
					  fftw_direction dir, int flags);
extern fftw_plan *fftwnd_new_plan_array(int rank);
extern fftw_plan *fftwnd_create_plans_generic(fftw_plan *plans,
					      int rank, const int *n,
					      fftw_direction dir, int flags);
extern fftw_plan *fftwnd_create_plans_specific(fftw_plan *plans,
					       int rank, const int *n,
					       const int *n_after,
					       fftw_direction dir, int flags,
					       fftw_complex *in, int istride,
					       fftw_complex *out, int ostride);
extern int fftwnd_work_size(int rank, const int *n, int flags, int ncopies);

extern void fftwnd_aux(fftwnd_plan p, int cur_dim,
		       fftw_complex *in, int istride,
		       fftw_complex *out, int ostride,
		       fftw_complex *work);
extern void fftwnd_aux_howmany(fftwnd_plan p, int cur_dim,
			       int howmany,
			       fftw_complex *in, int istride, int idist,
			       fftw_complex *out, int ostride, int odist,
			       fftw_complex *work);

/* wisdom prototypes */
enum fftw_wisdom_category {
     FFTW_WISDOM, RFFTW_WISDOM
};

extern int fftw_wisdom_lookup(int n, int flags, fftw_direction dir,
			      enum fftw_wisdom_category category,
			      int istride, int ostride,
			      enum fftw_node_type *type,
			      int *signature, int replace_p);
extern void fftw_wisdom_add(int n, int flags, fftw_direction dir,
			    enum fftw_wisdom_category cat,
			    int istride, int ostride,
			    enum fftw_node_type type,
			    int signature);

/* Private planner functions: */
extern double fftw_estimate_node(fftw_plan_node *p);
extern fftw_plan_node *fftw_make_node_notw(int size,
					const fftw_codelet_desc *config);
extern fftw_plan_node *fftw_make_node_real2hc(int size,
					const fftw_codelet_desc *config);
extern fftw_plan_node *fftw_make_node_hc2real(int size,
					const fftw_codelet_desc *config);
extern fftw_plan_node *fftw_make_node_twiddle(int n,
					 const fftw_codelet_desc *config,
					      fftw_plan_node *recurse,
					      int flags);
extern fftw_plan_node *fftw_make_node_hc2hc(int n,
					    fftw_direction dir,
					 const fftw_codelet_desc *config,
					    fftw_plan_node *recurse,
					    int flags);
extern fftw_plan_node *fftw_make_node_generic(int n, int size,
					      fftw_generic_codelet *codelet,
					      fftw_plan_node *recurse,
					      int flags);
extern fftw_plan_node *fftw_make_node_rgeneric(int n, int size,
					       fftw_direction dir,
					       fftw_rgeneric_codelet * codelet,
					       fftw_plan_node *recurse,
					       int flags);
extern int fftw_factor(int n);
extern fftw_plan_node *fftw_make_node(void);
extern fftw_plan fftw_make_plan(int n, fftw_direction dir,
				fftw_plan_node *root, int flags,
				enum fftw_node_type wisdom_type,
				int wisdom_signature);
extern void fftw_use_plan(fftw_plan p);
extern void fftw_use_node(fftw_plan_node *p);
extern void fftw_destroy_plan_internal(fftw_plan p);
extern fftw_plan fftw_pick_better(fftw_plan p1, fftw_plan p2);
extern fftw_plan fftw_lookup(fftw_plan *table, int n, int flags);
extern void fftw_insert(fftw_plan *table, fftw_plan this_plan, int n);
extern void fftw_make_empty_table(fftw_plan *table);
extern void fftw_destroy_table(fftw_plan *table);
extern void fftw_complete_twiddle(fftw_plan_node *p, int n);

extern fftw_plan_node *fftw_make_node_rader(int n, int size,
					    fftw_direction dir,
					    fftw_plan_node *recurse,
					    int flags);
extern fftw_rader_data *fftw_rader_top;

/****************************************************************************/
/*                           Floating Point Types                           */
/****************************************************************************/

/*
 * We use these definitions to make it easier for people to change
 * FFTW to use long double and similar types. You shouldn't have to
 * change this just to use float or double. 
 */

/*
 * Change this if your floating-point constants need to be expressed
 * in a special way.  For example, if fftw_real is long double, you
 * will need to append L to your fp constants to make them of the
 * same precision.  Do this by changing "x" below to "x##L". 
 */
#define FFTW_KONST(x) ((fftw_real) x)

#define FFTW_TRIG_SIN sin
#define FFTW_TRIG_COS cos
typedef double FFTW_TRIG_REAL;	/* the argument type for sin and cos */

#define FFTW_K2PI FFTW_KONST(6.2831853071795864769252867665590057683943388)

/****************************************************************************/
/*                               gcc/x86 hacks                              */
/****************************************************************************/

/*
 * gcc 2.[78].x and x86 specific hacks.  These macros align the stack
 * pointer so that the double precision temporary variables in the
 * codelets will be aligned to a multiple of 8 bytes (*way* faster on
 * pentium and pentiumpro)
 */
#ifdef __GNUC__
#ifdef __i386__
#ifdef FFTW_ENABLE_I386_HACKS
#ifndef FFTW_ENABLE_FLOAT
#define FFTW_USING_I386_HACKS
#define HACK_ALIGN_STACK_EVEN() {                        \
     if ((((long) (__builtin_alloca(0))) & 0x7)) __builtin_alloca(4);    \
}

#define HACK_ALIGN_STACK_ODD() {                         \
     if (!(((long) (__builtin_alloca(0))) & 0x7)) __builtin_alloca(4);   \
}

#ifdef FFTW_DEBUG_ALIGNMENT
#define ASSERT_ALIGNED_DOUBLE() {                        \
     double __foo;                                       \
     if ((((long) &__foo) & 0x7)) abort();                \
}
#endif

#endif
#endif
#endif
#endif

#ifndef HACK_ALIGN_STACK_EVEN
#define HACK_ALIGN_STACK_EVEN()
#endif
#ifndef HACK_ALIGN_STACK_ODD
#define HACK_ALIGN_STACK_ODD()
#endif
#ifndef ASSERT_ALIGNED_DOUBLE
#define ASSERT_ALIGNED_DOUBLE()
#endif

/****************************************************************************/
/*                                  Timers                                  */
/****************************************************************************/

/*
 * Here, you can use all the nice timers available in your machine.
 */

/*
 *
 Things you should define to include your own clock:
 
 fftw_time -- the data type used to store a time
 
 extern fftw_time fftw_get_time(void); 
 -- a function returning the current time.  (We have
 implemented this as a macro in most cases.)
 
 extern fftw_time fftw_time_diff(fftw_time t1, fftw_time t2);
 -- returns the time difference (t1 - t2).
 If t1 < t2, it may simply return zero (although this
 is not required).  (We have implemented this as a macro
 in most cases.)
 
 extern double fftw_time_to_sec(fftw_time t);
 -- returns the time t expressed in seconds, as a double.
 (Implemented as a macro in most cases.)
 
 FFTW_TIME_MIN -- a double-precision macro holding the minimum
 time interval (in seconds) for accurate time measurements.
 This should probably be at least 100 times the precision of
 your clock (we use even longer intervals, to be conservative).
 This will determine how long the planner takes to measure
 the speeds of different possible plans.
 
 Bracket all of your definitions with an appropriate #ifdef so that
 they will be enabled on your machine.  If you do add your own
 high-precision timer code, let us know (at fftw@theory.lcs.mit.edu).
 
 Only declarations should go in this file.  Any function definitions
 that you need should go into timer.c.
 */

/*
 * define a symbol so that we know that we have the fftw_time_diff
 * function/macro (it did not exist prior to FFTW 1.2) 
 */
#define FFTW_HAS_TIME_DIFF

/**********************************************
 *              SOLARIS
 **********************************************/
#if defined(HAVE_GETHRTIME)

/* we use the nanosecond virtual timer */
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

typedef hrtime_t fftw_time;

#define fftw_get_time() gethrtime()
#define fftw_time_diff(t1,t2) ((t1) - (t2))
#define fftw_time_to_sec(t) ((double) t / 1.0e9)

/*
 * a measurement is valid if it runs for at least
 * FFTW_TIME_MIN seconds.
 */
#define FFTW_TIME_MIN (1.0e-4)	/* for Solaris nanosecond timer */
#define FFTW_TIME_REPEAT 8

/**********************************************
 *        Pentium time stamp counter
 **********************************************/
#elif defined(__GNUC__) && defined(__i386__) && defined(FFTW_ENABLE_PENTIUM_TIMER)

/*
 * Use internal Pentium register (time stamp counter). Resolution
 * is 1/FFTW_CYCLES_PER_SEC seconds (e.g. 5 ns for Pentium 200 MHz).
 * (This code was contributed by Wolfgang Reimer)
 */

#ifndef FFTW_CYCLES_PER_SEC
#error "Must define FFTW_CYCLES_PER_SEC in fftw/config.h to use the Pentium cycle counter"
#endif

typedef unsigned long long fftw_time;

static __inline__ fftw_time read_tsc()
{
     struct {
	  long unsigned lo, hi;
     } counter;
     long unsigned sav_eax, sav_edx;
   __asm__("movl %%eax,%0":"=m"(sav_eax));
   __asm__("movl %%edx,%0":"=m"(sav_edx));
     __asm__("rdtsc");
   __asm__("movl %%eax,%0":"=m"(counter.lo));
   __asm__("movl %%edx,%0":"=m"(counter.hi));
   __asm__("movl %0,%%eax": : "m"(sav_eax):"eax");
   __asm__("movl %0,%%edx": : "m"(sav_edx):"edx");
     return *(fftw_time *) & counter;
}

#define fftw_get_time()  read_tsc()
#define fftw_time_diff(t1,t2) ((t1) - (t2))
#define fftw_time_to_sec(t) (((double) (t)) / FFTW_CYCLES_PER_SEC)
#define FFTW_TIME_MIN (1.0e-4)	/* for Pentium TSC register */

/************* generic systems having gettimeofday ************/
#elif defined(HAVE_GETTIMEOFDAY) || defined(HAVE_BSDGETTIMEOFDAY)
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#define FFTW_USE_GETTIMEOFDAY

typedef struct timeval fftw_time;

extern fftw_time fftw_gettimeofday_get_time(void);
extern fftw_time fftw_gettimeofday_time_diff(fftw_time t1, fftw_time t2);
#define fftw_get_time() fftw_gettimeofday_get_time()
#define fftw_time_diff(t1, t2) fftw_gettimeofday_time_diff(t1, t2)
#define fftw_time_to_sec(t) ((double)(t).tv_sec + (double)(t).tv_usec * 1.0E-6)

#ifndef FFTW_TIME_MIN
/* this should be fine on any system claiming a microsecond timer */
#define FFTW_TIME_MIN (1.0e-2)
#endif

/**********************************************
 *              MACINTOSH
 **********************************************/
#elif defined(HAVE_MAC_TIMER)

/*
 * By default, use the microsecond-timer in the Mac Time Manager.
 * Alternatively, by changing the following #if 1 to #if 0, you
 * can use the nanosecond timer available *only* on PCI PowerMacs. 
 */
#ifndef HAVE_MAC_PCI_TIMER	/* use time manager */

/*
 * Use Macintosh Time Manager routines (maximum resolution is about 20
 * microseconds). 
 */
typedef struct fftw_time_struct {
     unsigned long hi, lo;
} fftw_time;

extern fftw_time get_Mac_microseconds(void);

#define fftw_get_time() get_Mac_microseconds()

/* define as a function instead of a macro: */
extern fftw_time fftw_time_diff(fftw_time t1, fftw_time t2);

#define fftw_time_to_sec(t) ((t).lo * 1.0e-6 + 4294967295.0e-6 * (t).hi)

/* very conservative, since timer should be accurate to 20e-6: */
/* (although this seems not to be the case in practice) */
#define FFTW_TIME_MIN (5.0e-2)	/* for MacOS Time Manager timer */

#else				/* use nanosecond timer */

/* Use the nanosecond timer available on PCI PowerMacs. */

#include <DriverServices.h>

typedef AbsoluteTime fftw_time;
#define fftw_get_time() UpTime()
#define fftw_time_diff(t1,t2) SubAbsoluteFromAbsolute(t1,t2)
#define fftw_time_to_sec(t) (AbsoluteToNanoseconds(t).lo * 1.0e-9)

/* Extremely conservative minimum time: */
/* for MacOS PCI PowerMac nanosecond timer */
#define FFTW_TIME_MIN (5.0e-3)	

#endif				/* use nanosecond timer */

/**********************************************
 *              WINDOWS
 **********************************************/
#elif defined(HAVE_WIN32_TIMER)

#include <time.h>

typedef unsigned long fftw_time;
extern unsigned long GetPerfTime(void);
extern double GetPerfSec(double ticks);

#define fftw_get_time() GetPerfTime()
#define fftw_time_diff(t1,t2) ((t1) - (t2))
#define fftw_time_to_sec(t) GetPerfSec(t)

#define FFTW_TIME_MIN (5.0e-2)	/* for Win32 timer */

/**********************************************
 *              CRAY
 **********************************************/
#elif defined(_CRAYMPP)		/* Cray MPP system */

double SECONDR(void);		/* 
				 * I think you have to link with -lsci to
				 * get this 
				 */

typedef double fftw_time;
#define fftw_get_time() SECONDR()
#define fftw_time_diff(t1,t2) ((t1) - (t2))
#define fftw_time_to_sec(t) (t)

#define FFTW_TIME_MIN (1.0e-1)	/* for Cray MPP SECONDR timer */

/**********************************************
 *          VANILLA UNIX/ISO C SYSTEMS
 **********************************************/
/* last resort: use good old Unix clock() */
#else

#include <time.h>

typedef clock_t fftw_time;

#ifndef CLOCKS_PER_SEC
#ifdef sun
/* stupid sunos4 prototypes */
#define CLOCKS_PER_SEC 1000000
extern long clock(void);
#else				/* not sun, we don't know CLOCKS_PER_SEC */
#error Please define CLOCKS_PER_SEC
#endif
#endif

#define fftw_get_time() clock()
#define fftw_time_diff(t1,t2) ((t1) - (t2))
#define fftw_time_to_sec(t) (((double) (t)) / CLOCKS_PER_SEC)

/*
 * ***VERY*** conservative constant: this says that a
 * measurement must run for 200ms in order to be valid.
 * You had better check the manual of your machine
 * to discover if it can do better than this
 */
#define FFTW_TIME_MIN (2.0e-1)	/* for default clock() timer */

#endif				/* UNIX clock() */

/* take FFTW_TIME_REPEAT measurements... */
#ifndef FFTW_TIME_REPEAT
#define FFTW_TIME_REPEAT 4
#endif

/* but do not run for more than TIME_LIMIT seconds while measuring one FFT */
#ifndef FFTW_TIME_LIMIT
#define FFTW_TIME_LIMIT 2.0
#endif

#ifdef __cplusplus
}				/* extern "C" */

#endif				/* __cplusplus */

#endif				/* FFTW_INT_H */
