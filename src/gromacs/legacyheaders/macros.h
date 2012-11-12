/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _macros_h
#define _macros_h

#include "typedefs.h" /* for real definition only */

/* no extern "C" for this header because it only defines Macros */

/*
 * With the macros below you don't
 * have to use an index if you don't wan't to. You can eg. use
 * angle.C0[23] instead if angle.c[0][23].
 * In a similar fashion, you can use angle.AI[3] instead of
 * angle.a[0][3]
 */
#ifndef __cplusplus
#define AI  a[0]
#define AJ  a[1]
#define AK  a[2]
#define AL  a[3]
#define AM  a[4]
#define C0  c[0]
#define C1  c[1]
#define C2  c[2]

#ifndef min
#define min(a, b) (((a) < (b)) ? (a) : (b) )
#endif
#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b) )
#endif
#endif

/* This macro calculates the size of a array */
#define asize(a) (sizeof(a)/sizeof((a)[0]))

#endif  /* _macros_h */

