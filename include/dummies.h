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

#ifndef _dummies_h
#define _dummies_h

static char *SRCID_dummies_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"

typedef struct {
  int nprevdum; /* how many dummy particles are nonlocal */     
  int nnextdum;
  int *idxprevdum; /* index of nonlocal dummy particles */
  int *idxnextdum;
  int nprevconstr; /* how many constr. atoms are nonlocal */
  int nnextconstr;
  int *idxprevconstr; /* indices of nonlocal constructing atoms */
  int *idxnextconstr;
} t_comm_dummies;

extern void construct_dummies(FILE *log,rvec x[],t_nrnb *nrnb,
			      real dt,rvec v[],t_idef *idef,
			      t_graph *graph,t_commrec *cr,
			      matrix box,t_comm_dummies *dummycomm);
/* Create positions of dummy atoms based on surrounding atoms.
 */
 
extern void spread_dummy_f(FILE *log,rvec x[],rvec f[],
			   t_nrnb *nrnb,t_idef *idef,
			   t_comm_dummies *dummycomm,t_commrec *cr);
/* Spread the force operating on the dummy atoms on the surrounding atoms.
 */

#endif

