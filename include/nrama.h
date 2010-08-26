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

#ifndef _nrama_h
#define _nrama_h

#include "typedefs.h"
#include "statutil.h"
#include "mshift.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  gmx_bool bShow;
  char *label;
  int  iphi,ipsi; /* point in the dih array of xr... */
} t_phipsi;

typedef struct {
  atom_id ai[4];
  int     mult;
  real    phi0;
  real    ang;
} t_dih;

typedef struct {
  int       ndih;
  t_dih     *dih;
  int       npp;
  t_phipsi  *pp;
  t_trxstatus *traj;
  int       natoms;
  int       amin,amax;
  real      t;
  rvec      *x;
  matrix    box;
  t_idef    *idef;
  int       ePBC;
  output_env_t oenv;
} t_xrama;

t_topology *init_rama(const output_env_t oenv, const char *infile,
                             const char *topfile, t_xrama *xr,int mult);

gmx_bool new_data(t_xrama *xr);

#ifdef __cplusplus
}
#endif

#endif	/* _nrama_h */
