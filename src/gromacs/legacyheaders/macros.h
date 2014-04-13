/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef _macros_h
#define _macros_h

#include <stdlib.h>

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
#define asize(a) ((int)(sizeof(a)/sizeof((a)[0])))

#endif  /* _macros_h */
