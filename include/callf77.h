/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
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
