/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Getting the Right Output Means no Artefacts in Calculating Stuff
 */

#ifndef _rbin_h
#define _rbin_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sysstuff.h"
#include "typedefs.h"
#include "network.h"
	
typedef struct {
  int    nreal;
  int    maxreal;
  double *rbuf;
} t_bin;

extern t_bin *mk_bin(void);
/* Create a real bin */

extern void reset_bin(t_bin *b);
/* Reset number of entries to zero */

extern int add_binr(FILE *log,t_bin *b,int nr,real r[]);
extern int add_bind(FILE *log,t_bin *b,int nr,double r[]);
/* Add reals to the bin. Returns index */

extern void sum_bin(t_bin *b,t_commrec *cr);
/* Globally sum the reals in the bin */

extern void extract_binr(t_bin *b,int index,int nr,real r[]);
extern void extract_bind(t_bin *b,int index,int nr,double r[]);
/* Extract values from the bin, starting from index (see add_bin) */

#endif	/* _rbin_h */
