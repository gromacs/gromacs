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

extern void reset_index(t_mat *m);

extern void swap_rows(t_mat *m,int isw,int jsw);

extern void set_mat_entry(t_mat *m,int i,int j,real val);

extern void done_mat(t_mat **m);

extern real row_energy(int n1,int row,real *mat,int m_ind[]);

extern real mat_energy(t_mat *mat);

extern void swap_mat(t_mat *m);

extern void rms_dist(char *fn,t_mat *m);

extern t_clustid *new_clustid(int n1);

#endif



