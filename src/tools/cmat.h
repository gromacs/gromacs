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
 * Great Red Oystrich Makes All Chemists Sane
 */

#ifndef _cmat_h
#define _cmat_h

static char *SRCID_cmat_h = "$Id$";

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
  bool b1D;
  real emat,maxrms,sumrms;
  real *erow;
  real **mat;
} t_mat;

/* The matrix is indexed using the matrix index */
#define EROW(m,i)  m->erow[i]

extern real **mk_matrix(int n1,bool b1D);

extern void done_matrix(int n1,real ***m);

extern t_mat *init_mat(int n1,bool b1D);

extern void enlarge_mat(t_mat *m,int deltan);

extern void reset_index(t_mat *m);

extern void swap_rows(t_mat *m,int isw,int jsw);

extern void set_mat_entry(t_mat *m,int i,int j,real val);

extern void done_mat(t_mat **m);

extern real row_energy(int n1,int row,real *mat,int m_ind[]);

extern real mat_energy(t_mat *mat);

extern void swap_mat(t_mat *m);

extern void low_rmsd_dist(char *fn,real maxrms,int nn,real **mat);

extern void rmsd_distribution(char *fn,t_mat *m);

extern t_clustid *new_clustid(int n1);

#endif



