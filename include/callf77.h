/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Giving Russians Opium May Alter Current Situation
 */

#ifndef _callf77_h
#define _callf77_h

static char *SRCID_callf77_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef USE_FORTRAN

#define FUNC(name,NAME) F77_FUNC(name,NAME)
#define SCALARG(name) *name
#define SCAL(name) &(name)

/* define f77 name mangling - we dont need any
 * special macros for names with underscores since
 * no such identifiers exist in gromacs right now.
 * If you add one you should include the definition
 * of F77_NAME_EXTRA_UNDERSCORE below and create
 * the macro F77_FUNC_(name,NAME).
 */
#ifdef F77_NAME_LOWERCASE
#  define F77_FUNC(name,NAME)     name
#  ifdef F77_NAME_EXTRA_UNDERSCORE
#    define F77_FUNC_(name,NAME)  name ## _
#  else
#    define F77_FUNC_(name,NAME)  name
#  endif
#elif defined F77_NAME_LOWERCASE_UNDERSCORE
#  define F77_FUNC(name,NAME)     name ## _
#  ifdef F77_NAME_EXTRA_UNDERSCORE
#    define F77_FUNC_(name,NAME)  name ## __
#  else
#    define F77_FUNC_(name,NAME)  name ## _
#  endif
#elif defined F77_NAME_UPPERCASE
#  define F77_FUNC(name,NAME)     NAME
#  ifdef F77_NAME_EXTRA_UNDERSCORE
#    define F77_FUNC_(name,NAME)  NAME ## _
#  else
#    define F77_FUNC_(name,NAME)  NAME
#  endif
#elif defined F77_NAME_UPPERCASE_UNDERSCORE
#  define F77_FUNC(name,NAME)     NAME ## _
#  ifdef F77_NAME_EXTRA_UNDERSCORE
#    define F77_FUNC_(name,NAME)  NAME ## __
#  else
#    define F77_FUNC_(name,NAME)  NAME ## _
#  endif
#endif

#else /* Use C */

#define FUNC(name,NAME) name
#define SCALARG(name) name
#define SCAL(name) name

#endif

#endif
