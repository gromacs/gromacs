/* -*- C -*- */
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

/* fftw.h -- system-wide definitions */
/* $Id$ */

/* configuration options (guessed by configure) */

/* Define to empty if the keyword does not work.  */
/* #undef const */

/* Define if you have the gettimeofday function.  */
/* #undef HAVE_GETTIMEOFDAY */

/* Define if you have the BSDgettimeofday function.  */
/* #undef HAVE_BSDGETTIMEOFDAY */

/* Define if you have the <sys/time.h> header file.  */
/* #undef HAVE_SYS_TIME_H */

/* Define if you have the <unistd.h> header file.  */
/* #undef HAVE_UNISTD_H */

/* Define if you have the <getopt.h> header file.  */
/* #undef HAVE_GETOPT_H */

/* Define if you have the <malloc.h> header file */
/* #undef HAVE_MALLOC_H */

/* Define if you have gethrtime() a la Solaris 2 */
/* #undef HAVE_GETHRTIME */

/* Define if you have getopt() */
/* #undef HAVE_GETOPT */

/* Define if you have isnan() */
/* #undef HAVE_ISNAN */

/* Define for enabling the high resolution Pentium timer */
/* #undef FFTW_ENABLE_PENTIUM_TIMER */

/*
 * When using FFTW_ENABLE_PENTIUM_TIMER, set FFTW_CYCLES_PER_SEC 
 * to your real CPU clock speed! 
 */
/* This is for 200 MHz */
/* #define FFTW_CYCLES_PER_SEC 200000000L */

/*
 * Define to enable a gcc/x86 specific hack that aligns
 * the stack to an 8-byte boundary 
 */
/* #undef FFTW_ENABLE_I386_HACKS */

/* Define to enable extra runtime checks for debugging. */
/* #undef FFTW_DEBUG */

/*
 * Define to enable extra runtime checks for the alignment of variables
 * in the codelets (causes coredump for misaligned double on x86). 
 */
/* #undef FFTW_DEBUG_ALIGNMENT */

#define FFTW_VERSION "2.0.1" 

/* Use Win32 high-resolution timer */
#if defined(__WIN32__) || defined(WIN32) || defined(_WINDOWS)
#define HAVE_WIN32_TIMER
#endif

/* Use MacOS Time Manager timer */
#if defined(MAC) || defined(macintosh)
#define HAVE_MAC_TIMER

/* Define to use nanosecond timer on PCI PowerMacs: */
/* #undef HAVE_MAC_PCI_TIMER */
#endif
