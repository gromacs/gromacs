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
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _sortwater_h
#define _sortwater_h

static char *SRCID_sortwater_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"

extern void randwater(int astart,int nwater,int nwatom,
		      rvec x[],rvec v[],int *seed);
/* Randomize the order of nwater molecules of length nwatom, the
 * first atom of which is at astart.
 * If v is not NULL it will be shuffled along
 */


extern void sortwater(int astart,int nwater,int nwatom,rvec x[],rvec v[]);
/* Sort the order of nwater molecules of length nwatom on X coordinate
 * If v is not NULL it will be shuffled along
 */

extern void mkcompact(int astart,int nwater,int nwatom,rvec x[],rvec v[],
		      int nnode,matrix box);
/* Make compact subboxes */

#endif
