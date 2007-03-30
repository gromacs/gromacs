/*
 * $Id$
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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef _constr_h
#define _constr_h

/* LINCS stuff */
typedef struct {
  int  nc;       /* the number of constraints */
  int  nc_alloc; /* the number we allocated memory for */
  int  nflexcon; /* the number of flexible constraints */
  int  ncc;      /* the number of constraint connections */
  int  ncc_alloc;/* the number we allocated memory for */
  real matlam;   /* the FE lambda value used for filling blc and blcc */
  real *bllen0;  /* the reference distance in topology A */
  real *ddist;   /* the reference distance in top B - the r.d. in top A */
  int  *bla;     /* the atom pairs involved in the constraints */
  real *blc;     /* 1/sqrt(invmass1 + invmass2) */
  int  *blnr;    /* index into blbnb and blcc */
  int  *blbnb;   /* list of bond connections */
  real *blcc;    /* bond coupling coefficient matrix */
  real *bllen;   /* the reference bond length */
  /* arrays for temporary storage in the LINCS algorithm */
  rvec *tmpv;
  real *tmpncc;
  real *tmp1;
  real *tmp2;
  real *tmp3;
  real *lambda;  /* the Lagrange multipliers */
} t_lincsdata;

/* All the constraint data */
typedef struct {
  int         nflexcon;     /* The number of flexible constraints */
  t_lincsdata *lincsd;      /* LINCS data                         */
  int         nblocks;      /* The number of SHAKE blocks         */
  int         *sblock;      /* The SHAKE blocks                   */
  int         maxwarn;      /* The maximum number of warnings     */
  int         warncount_lincs;
  int         warncount_settle;
} gmx_constr_t;

#endif
