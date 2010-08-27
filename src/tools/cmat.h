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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _cmat_h
#define _cmat_h

#include "typedefs.h"

typedef struct {
  int  i,j;
  real dist;
} t_dist;

typedef struct {
  int  conf,clust;
} t_clustid;

typedef struct {
  int  n1,nn;
  int  *m_ind;
  gmx_bool b1D;
  real emat,minrms,maxrms,sumrms;
  real *erow;
  real **mat;
} t_mat;

/* The matrix is indexed using the matrix index */
#define EROW(m,i)  m->erow[i]

extern t_mat *init_mat(int n1,gmx_bool b1D);

extern void enlarge_mat(t_mat *m,int deltan);

extern void reset_index(t_mat *m);

extern void swap_rows(t_mat *m,int isw,int jsw);

extern void set_mat_entry(t_mat *m,int i,int j,real val);

extern void done_mat(t_mat **m);

extern real row_energy(int n1,int row,real *mat);

extern real mat_energy(t_mat *mat);

extern void swap_mat(t_mat *m);

extern void low_rmsd_dist(const char *fn,real maxrms,int nn,real **mat,
                          const output_env_t oenv);

extern void rmsd_distribution(const char *fn,t_mat *m,const output_env_t oenv);

extern t_clustid *new_clustid(int n1);

#endif
