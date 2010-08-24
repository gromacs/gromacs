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

#ifndef _sortwater_h
#define _sortwater_h

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void randwater(int astart,int nwater,int nwatom,
		      rvec x[],rvec v[],int *seed);
/* Randomize the order of nwater molecules of length nwatom, the
 * first atom of which is at astart.
 * If v is not NULL it will be shuffled along
 * IS NOT THREAD SAFE 
 */


void sortwater(int astart,int nwater,int nwatom,rvec x[],rvec v[]);
/* Sort the order of nwater molecules of length nwatom on X coordinate
 * If v is not NULL it will be shuffled along
 * IS NOT THREAD SAFE 
 */

void mkcompact(int astart,int nwater,int nwatom,rvec x[],rvec v[],
		      int nnode,matrix box);
/* Make compact subboxes 
 * IS NOT THREAD SAFE  */

#ifdef __cplusplus
}
#endif

#endif
