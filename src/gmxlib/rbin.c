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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "main.h"
#include "network.h"
#include "rbin.h"
#include "smalloc.h"
	
t_bin *mk_bin(void)
{
  t_bin *b;
  
  snew(b,1);
  
  return b;
}

void destroy_bin(t_bin *b)
{
  if (b->maxreal > 0) {
    sfree(b->rbuf);
  }
  
  sfree(b);
}

void reset_bin(t_bin *b)
{
  b->nreal = 0;
}

int add_binr(t_bin *b,int nr,real r[])
{
#define MULT 4
  int    i,rest,index;
  double *rbuf;
  
  if (b->nreal+nr > b->maxreal) {
    b->maxreal=b->nreal+nr;
    rest=b->maxreal % MULT;
    if (rest != 0)
      b->maxreal+=MULT-rest;
    srenew(b->rbuf,b->maxreal);
  }
  /* Copy pointer */
  rbuf=b->rbuf+b->nreal;
  for(i=0; (i<nr); i++)
    rbuf[i]=r[i];
    
  index=b->nreal;
  b->nreal+=nr;
  
  return index;
}

int add_bind(t_bin *b,int nr,double r[])
{
#define MULT 4
  int    i,rest,index;
  double *rbuf;
  
  if (b->nreal+nr > b->maxreal) {
    b->maxreal=b->nreal+nr;
    rest=b->maxreal % MULT;
    if (rest != 0)
      b->maxreal+=MULT-rest;
    srenew(b->rbuf,b->maxreal);
  }
  /* Copy pointer */
  rbuf=b->rbuf+b->nreal;
  for(i=0; (i<nr); i++)
    rbuf[i]=r[i];
    
  index=b->nreal;
  b->nreal+=nr;
  
  return index;
}

void sum_bin(t_bin *b,t_commrec *cr)
{
  int i;
  
  for(i=b->nreal; (i<b->maxreal); i++)
    b->rbuf[i]=0;
  gmx_sumd(b->maxreal,b->rbuf,cr);
}

void extract_binr(t_bin *b,int index,int nr,real r[])
{
  int    i;
  double *rbuf;
  
  rbuf = b->rbuf+index;
  for(i=0; (i<nr); i++)
    r[i]=rbuf[i];
}

void extract_bind(t_bin *b,int index,int nr,double r[])
{
  int    i;
  double *rbuf;
  
  rbuf = b->rbuf+index;
  for(i=0; (i<nr); i++)
    r[i]=rbuf[i];
}

#ifdef DEBUGRBIN
int main(int argc,char *argv[])
{
  t_commrec *cr;
  t_bin     *rb;
  double    *r;
  rvec      *v;
  int       k,i,ni,mi,n,m;

  cr=init_par(&argc,argv);
  n=strtol(argv[1],NULL,10);
  m=strtol(argv[2],NULL,10);
  fprintf(stdlog,"n=%d\n",n);
  rb=mk_bin();
  snew(r,n);
  snew(v,m);
  
  for(k=0; (k < 3); k++) {
    fprintf(stdlog,"\nk=%d\n",k);
    reset_bin(rb);
    
    for(i=0; (i<n); i++)
      r[i]=i+k;
    for(i=0; (i<m); i++) {
      v[i][XX]=4*i+k;
      v[i][YY]=4*i+k+1;
      v[i][ZZ]=4*i+k+2;
    }

    ni=add_bind(stdlog,rb,n,r);
    mi=add_binr(stdlog,rb,DIM*m,v[0]);
    
    sum_bin(rb,cr);
    
    extract_bind(rb,ni,n,r);
    extract_binr(rb,mi,DIM*m,v[0]);
  
    for(i=0; (i<n); i++)
      fprintf(stdlog,"r[%d] = %e\n",i,r[i]);
    for(i=0; (i<m); i++)
      fprintf(stdlog,"v[%d] = (%e,%e,%e)\n",i,v[i][XX],v[i][YY],v[i][ZZ]);
  }
  fflush(stdlog);
  
  return 0;
}
#endif
