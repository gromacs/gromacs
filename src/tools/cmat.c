#include <math.h>
#include "cmat.h"
#include "smalloc.h"
#include "macros.h"
#include "xvgr.h"

real **mk_matrix(int n1,bool b1D)
{
  real **m;
  int  i;
  
  snew(m,n1);
  if (b1D)
    snew(m[0],n1*n1); 
  
  for(i=0; (i<n1); i++) {
    if (b1D)
      m[i] = &(m[0][i*n1]);
    else
      snew(m[i],n1);
  }
  return m;
}

void done_matrix(int n1,real ***m)
{
  int i;
  
  for(i=0; (i<n1); i++)
    sfree((*m)[i]);
  sfree(*m);
  *m = NULL;
}

t_mat *init_mat(int n1,bool b1D)
{
  t_mat *m;
  int   i;
  
  snew(m,1);
  m->n1   = n1;
  m->nn   = 0;
  m->b1D  = b1D;
  m->emat = 0;
  m->maxrms = 0;
  m->sumrms = 0;
  m->nn   = 0;
  m->mat  = mk_matrix(n1,b1D);
  
  snew(m->erow,n1);
  snew(m->m_ind,n1);
  reset_index(m);
    
  return m;
}

void reset_index(t_mat *m)
{
  int i;
  
  for(i=0; (i<m->n1); i++)
    m->m_ind[i] = i;
}

void set_mat_entry(t_mat *m,int i,int j,real val)
{
  m->mat[i][j] = m->mat[j][i] = val;
  m->maxrms    = max(m->maxrms,val);
  m->sumrms   += val;
  m->nn        = max(m->nn,max(j+1,i+1));
}
  
void done_mat(t_mat **m)
{
  int i;
  
  done_matrix((*m)->n1,&((*m)->mat));  
  sfree((*m)->m_ind);
  sfree((*m)->erow);
  sfree(*m);
  *m = NULL;
}

real row_energy(int n1,int row,real *mat,int m_ind[])
{
  real re = 0;
  real n_1;
  int  i,ii;
  
  n_1 = 1.0/n1;
  for(i=0; (i<n1); i++) {
    ii  = m_ind[i];
    re += (abs(ii-row)*n_1)*mat[ii];
  }
  return re;
}

real mat_energy(t_mat *m)
{
  real re,retot;
  int  j,jj;
  
  retot = 0;
  for(j=0; (j<m->nn); j++) {
    jj = m->m_ind[j];
    re = row_energy(m->nn,j,m->mat[jj],m->m_ind);
    m->erow[j] = re;
    retot += re;
  }
  m->emat = retot/m->nn;
  return m->emat;
}

void swap_rows(t_mat *m,int isw,int jsw)
{
  real *tmp,ttt;
  int  i;

  /* Swap rows */
  tmp         = m->mat[isw];
  m->mat[isw] = m->mat[jsw];
  m->mat[jsw] = tmp;
  /* Swap columns */
  for(i=0; (i<m->nn); i++) {
    ttt            = m->mat[isw][i];
    m->mat[isw][i] = m->mat[jsw][i];
    m->mat[jsw][i] = ttt;
  }
}


void swap_mat(t_mat *m)
{
  t_mat *tmp;
  int   i,j;
  
  tmp = init_mat(m->nn,FALSE);
  for(i=0; (i<m->nn); i++)
    for(j=0; (j<m->nn); j++)
      tmp->mat[m->m_ind[i]][m->m_ind[j]] = m->mat[i][j]; 
  /*tmp->mat[i][j] =  m->mat[m->m_ind[i]][m->m_ind[j]]; */
  for(i=0; (i<m->nn); i++)
    for(j=0; (j<m->nn); j++)
      m->mat[i][j] = tmp->mat[i][j];
  done_mat(&tmp);
}

void rms_dist(char *fn,t_mat *rms)
{
  FILE   *fp;
  int    i,j,*histo;
  real   fac;
  
  fac = 100/rms->maxrms;
  snew(histo,101);
  for(i=0; (i<rms->nn); i++) 
    for(j=i+1; (j<rms->nn); j++)
      histo[(int)(fac*rms->mat[i][j])]++;
      
  fp = xvgropen(fn,"RMS Distribution","RMS (nm)","a.u.");
  for(i=0; (i<101); i++)
    fprintf(fp,"%10g  %10d\n",i/fac,histo[i]);
  fclose(fp);
  sfree(histo);
  xvgr_file(fn,NULL);
}

t_clustid *new_clustid(int n1)
{
  t_clustid *c;
  int       i;
  
  snew(c,n1);
  for(i=0; (i<n1); i++) {
    c[i].conf = i;
    c[i].clust = i;
  }
  return c;  
}

