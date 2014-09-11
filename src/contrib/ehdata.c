/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <math.h>
#include "typedefs.h"
#include "gromacs/utility/smalloc.h"
#include "macros.h"
#include "gromacs/utility/fatalerror.h"
#include "random.h"
#include "gromacs/fileio/strdb.h"
#include "gromacs/utility/futil.h"
#include "ehdata.h"

typedef struct {
  int  nener,nomega,nq;     /* Number of entries in the table      */
  real *ener;               /* Energy values                       */
  real **omega;             /* Omega values (energy dependent)     */
  real ***prob,***q;        /* Probability and Q                   */
} t_pq_inel;

typedef struct {
  int  nener,n2Ddata;       /* Number of entries in the table      */
  real *ener;               /* Energy values                       */
  real **prob,**data;       /* Probability and data (energy loss)  */
} t_p2Ddata;

static int realcomp(const void *a,const void *b)
{
  real *ra,*rb;
  
  ra = (real *)a;
  rb = (real *)b;
  if (ra < rb)
    return -1;
  else if (ra > rb)
    return 1;
  else 
    return 0;
}

static t_p2Ddata *read_p2Ddata(char *fn)
{
  FILE    *fp;
  t_p2Ddata *p2Ddata;
  int     i,j;
  double  e,p,o;
  
  fprintf(stdout,"Going to read %s\n",fn);
  fp = gmx_ffopen(fn,"r");

  /* Allocate memory and set constants */
  snew(p2Ddata,1);
  if (fscanf(fp,"%d%d",&p2Ddata->nener,&p2Ddata->n2Ddata) != 2)
    gmx_fatal(FARGS,"I need two integers: nener, n in file %s",fn);
  
  snew(p2Ddata->ener,p2Ddata->nener);
  snew(p2Ddata->prob,p2Ddata->nener);
  snew(p2Ddata->data,p2Ddata->nener);
  
  /* Double loop to read data */
  for(i=0; (i<p2Ddata->nener); i++) {
    fprintf(stderr,"\rEnergy %d/%d",i+1,p2Ddata->nener);
    snew(p2Ddata->prob[i],p2Ddata->n2Ddata);
    snew(p2Ddata->data[i],p2Ddata->n2Ddata);
    
    for(j=0; (j<p2Ddata->n2Ddata); j++) {
      fscanf(fp,"%lf%lf%lf",&e,&p,&o);

      /* Consistency check */
      if (j==0)
	p2Ddata->ener[i] = e;
      else if (fabs(p2Ddata->ener[i]-e) > 1e-6*e)
	gmx_fatal(FARGS,"Inconsistent energy %f i=%d j=%d",e,i,j);
      p2Ddata->prob[i][j] = p;
      p2Ddata->data[i][j] = o;
    }
    /* There is some noise on the data, take it away by sorting,
     * because otherwise binary search does not work.
     * This is equivalent to shifting in the data slightly along the X-axis
     * but better than linear search with the "real" data.
     */
    qsort(p2Ddata->data[i],p2Ddata->n2Ddata,sizeof(p2Ddata->data[0][0]),
	  realcomp);
  }
  fprintf(stderr,"\n");
  
  gmx_ffclose(fp);
  
  return p2Ddata;
}

static t_pq_inel *read_pq(char *fn)
{
  FILE      *fp;
  t_pq_inel *pq;
  int       i,j,k;
  double    e,p,o,t;
  
  fprintf(stdout,"Going to read %s\n",fn);
  fp = gmx_ffopen(fn,"r");

  /* Allocate memory and set constants */
  snew(pq,1);
  if (fscanf(fp,"%d%d%d",&pq->nener,&pq->nomega,&pq->nq) != 3)
    gmx_fatal(FARGS,"I need three integers: nener, nomega, nq in file %s",fn);
  
  snew(pq->ener,pq->nener);
  snew(pq->omega,pq->nener);
  snew(pq->prob,pq->nener);
  snew(pq->q,pq->nener);
  
  /* Triple loop to read data */
  for(i=0; (i<pq->nener); i++) {
    fprintf(stderr,"\rEnergy %d/%d",i+1,pq->nener);
    snew(pq->prob[i],pq->nomega);
    snew(pq->q[i],pq->nomega);
    snew(pq->omega[i],pq->nomega);
    
    for(j=0; (j<pq->nomega); j++) {
      snew(pq->prob[i][j],pq->nq);
      snew(pq->q[i][j],pq->nq);
      
      for(k=0; (k<pq->nq); k++) {
	fscanf(fp,"%lf%lf%lf%lf",&e,&o,&p,&t);
	
	/* Consistency check */
	if ((j == 0) && (k == 0)) 
	  pq->ener[i] = e;
	else if (fabs(pq->ener[i]-e) > 1e-6*e)
	  gmx_fatal(FARGS,"Inconsistent energy %f i=%d j=%d k=%d",e,i,j,k);
	
	if (k == 0)
	  pq->omega[i][j] = o;
	else if (fabs(pq->omega[i][j]-o) > 1e-6*o)
	  gmx_fatal(FARGS,"Inconsistent omega %f i=%d j=%d k=%d",o,i,j,k);
	
	pq->prob[i][j][k] = p;
	pq->q[i][j][k] = t;
      }
    }
  }
  fprintf(stderr,"\n");
  
  gmx_ffclose(fp);
  
  return pq;
}

static int my_bsearch(real val,int ndata,real data[],gmx_bool bLower)
{
  int ilo,ihi,imed;

  if (val < data[0])
    return -1;
  ilo = 0; 
  ihi = ndata-1;
  while ((ihi - ilo) > 1) {
    imed = (ilo+ihi)/2;
    if (data[imed] > val) 
      ihi = imed;
    else
      ilo = imed;
  }
  /* Now val should be in between data[ilo] and data[ihi] */
  /* Decide which one is closest */
  if (bLower || ((val-data[ilo]) <= (data[ihi]-val)))
    return ilo;
  else
    return ihi;
}

static real interpolate2D(int nx,int ny,real *dx,real **data,
			  real x0,real fy,int nx0,int ny0)
{
  real fx;
  
  fx  = (x0-dx[nx0])/(dx[nx0+1]-dx[nx0]);
  
  return (fx*fy*data[nx0][ny0] + fx*(1-fy)*data[nx0][ny0+1] +
	  (1-fx)*fy*data[nx0+1][ny0] + (1-fx)*(1-fy)*data[nx0+1][ny0+1]); 
}

real get_omega(real ekin,int *seed,FILE *fp,char *fn)
{
  static t_p2Ddata *p2Ddata = NULL;
  real r,ome,fx,fy;
  int  eindex,oindex;
  
  if (p2Ddata == NULL) 
    p2Ddata = read_p2Ddata(fn);
  
  /* Get energy index by binary search */
  if ((eindex = my_bsearch(ekin,p2Ddata->nener,p2Ddata->ener,TRUE)) >= 0) {
#ifdef DEBUG
    if (eindex >= p2Ddata->nener)
      gmx_fatal(FARGS,"eindex (%d) out of range (max %d) in get_omega",
		  eindex,p2Ddata->nener);
#endif

    /* Start with random number */
    r = rando(seed);
    
    /* Do binary search in the energy table */
    if ((oindex = my_bsearch(r,p2Ddata->n2Ddata,p2Ddata->prob[eindex],TRUE)) >= 0) {
#ifdef DEBUG
      if (oindex >= p2Ddata->n2Ddata)
	gmx_fatal(FARGS,"oindex (%d) out of range (max %d) in get_omega",
		    oindex,p2Ddata->n2Ddata);
#endif

      fy = ((r-p2Ddata->prob[eindex][oindex])/
	    (p2Ddata->prob[eindex][oindex+1]-p2Ddata->prob[eindex][oindex]));
      ome = interpolate2D(p2Ddata->nener,p2Ddata->n2Ddata,p2Ddata->ener,
			  p2Ddata->data,ekin,fy,
			  eindex,oindex);
      /* ome = p2Ddata->data[eindex][oindex];*/
      
      if (fp) 
	fprintf(fp,"%8.3f  %8.3f\n",ome,r);
      
      return ome;
    }
  }
  return 0;
}

real get_theta_el(real ekin,int *seed,FILE *fp,char *fn)
{
  static t_p2Ddata *p2Ddata = NULL;
  real r,theta;
  int  eindex,tindex;
  
  if (p2Ddata == NULL) 
    p2Ddata = read_p2Ddata(fn);
  
  /* Get energy index by binary search */
  if ((eindex = my_bsearch(ekin,p2Ddata->nener,p2Ddata->ener,TRUE)) >= 0) {
  
    /* Start with random number */
    r = rando(seed);
    
    /* Do binary search in the energy table */
    if ((tindex = my_bsearch(r,p2Ddata->n2Ddata,p2Ddata->prob[eindex],FALSE)) >= 0) {
  
      theta = p2Ddata->data[eindex][tindex];
      
      if (fp) 
	fprintf(fp,"%8.3f  %8.3f\n",theta,r);
      
      return theta;
    }
  }
  return 0;
}

real get_q_inel(real ekin,real omega,int *seed,FILE *fp,char *fn)
{
  static t_pq_inel *pq = NULL;
  int    eindex,oindex,tindex;
  real   r,theta;
  
  if (pq == NULL)
    pq = read_pq(fn);

  /* Get energy index by binary search */
  if ((eindex = my_bsearch(ekin,pq->nener,pq->ener,TRUE)) >= 0) {
#ifdef DEBUG
    if (eindex >= pq->nener)
      gmx_fatal(FARGS,"eindex out of range (%d >= %d)",eindex,pq->nener);
#endif
      
    /* Do binary search in the energy table */
    if ((oindex = my_bsearch(omega,pq->nomega,pq->omega[eindex],FALSE)) >= 0) {
#ifdef DEBUG
      if (oindex >= pq->nomega)
	gmx_fatal(FARGS,"oindex out of range (%d >= %d)",oindex,pq->nomega);
#endif
      
      /* Start with random number */
      r = rando(seed);
      
      if ((tindex = my_bsearch(r,pq->nq,pq->prob[eindex][oindex],FALSE)) >= 0) {
#ifdef DEBUG
	if (tindex >= pq->nq)
	  gmx_fatal(FARGS,"tindex out of range (%d >= %d)",tindex,pq->nq);
#endif
	
	theta = pq->q[eindex][oindex][tindex];
  
	if (fp)
	  fprintf(fp,"get_q_inel: %8.3f  %8.3f\n",theta,r);
	
	return theta;
      }
    }
  }
  return 0;
}

static int read_cross(char *fn,real **ener,real **cross,real factor)
{
  char   **data=NULL;
  double e,c;
  int    i,j,n;
  
  fprintf(stdout,"Going to read %s\n",fn);
  n = get_file(fn,&data);

  /* Allocate memory */
  snew(*cross,n);
  snew(*ener,n);
  for(i=j=0; (i<n); i++) {
    if (sscanf(data[i],"%lf%lf",&e,&c) == 2) {
      (*ener)[j] = e;
      (*cross)[j] = c*factor;
      j++;
    }
    sfree(data[i]);
  }
  sfree(data);
  if (j != n)
    fprintf(stderr,"There were %d empty lines in file %s\n",n-j,fn);
  
  return j;
}

real cross_inel(real ekin,real rho,char *fn)
{
  static real *ener  = NULL;
  static real *cross = NULL;
  static int  ninel;
  int eindex;
  
  /* Read data at first call, convert A^2 to nm^2 */
  if (cross == NULL)
    ninel = read_cross(fn,&ener,&cross,rho*0.01);
  
  /* Compute index with binary search */
  if ((eindex = my_bsearch(ekin,ninel,ener,FALSE)) >= 0) {
#ifdef DEBUG
    if (eindex >= ninel)
      gmx_fatal(FARGS,"ekin = %f, ener[0] = %f, ener[%d] = %f",ekin,
		  ener[0],ener[ninel-1]);
#endif
    return cross[eindex];
  }
  
  return 0;
}

real cross_el(real ekin,real rho,char *fn)
{
  static real *ener  = NULL;
  static real *cross = NULL;
  static int  nel;
  int eindex;
  
  /* Read data at first call, convert A^2 to nm^2  */
  if (cross == NULL) 
    nel = read_cross(fn,&ener,&cross,rho*0.01);
  
  /* Compute index with binary search */
  if ((eindex = my_bsearch(ekin,nel,ener,FALSE)) >= 0) {
#ifdef DEBUG
    if (eindex >= nel)
      gmx_fatal(FARGS,"ekin = %f, ener[0] = %f, ener[%d] = %f",ekin,
		  ener[0],ener[nel-1]);
#endif

    return cross[eindex];
  }
  return 0;
}

real band_ener(int *seed,FILE *fp,char *fn)
{
  static real *ener  = NULL;
  static real *prob  = NULL;
  static int  nener;
  int  eindex;
  real r;
  
  /* Read data at first call, misuse read_cross function */
  if (prob == NULL)
    nener = read_cross(fn,&ener,&prob,1.0);
  
  r = rando(seed);
  
  if ((eindex = my_bsearch(r,nener,prob,FALSE)) >= 0) {
#ifdef DEBUG
    if (eindex >= nener)
      gmx_fatal(FARGS,"r = %f, prob[0] = %f, prob[%d] = %f",r,
		  prob[0],prob[nener-1]);
#endif

    if (fp)
      fprintf(fp,"%8.3f  %8.3f\n",ener[eindex],r);
    
    return ener[eindex];
  }
  return 0;
}

static void test_omega(FILE *fp,int *seed)
{
  int i;

  fprintf(fp,"Testing the energy loss tables\n");
  for(i=0; (i<1000); i++) {
    (void) get_omega(500*rando(seed),seed,fp,NULL);
  }
}

static void test_q_inel(FILE *fp,int *seed)
{
  int i;
  
  fprintf(fp,"Testing the energy/omega dependent inelastic scattering q tables\n");
  for(i=0; (i<1000); i++) {
    (void) get_q_inel(500*rando(seed),400*rando(seed),seed,fp,NULL);
  }
}

static void test_theta_el(FILE *fp,int *seed)
{
  int i;
  
  fprintf(fp,"Testing the energy dependent elastic scattering theta tables\n");
  for(i=0; (i<1000); i++) {
    (void) get_theta_el(500*rando(seed),seed,fp,NULL);
  }
}

static void test_sigma_el(FILE *fp,real rho)
{
  int  i;

  fprintf(fp,"Testing the elastic cross sections table\n");
  for(i=0; (i<500); i++) {
    fprintf(fp,"%3d  %8.3f\n",i,cross_el(i,rho,NULL));
  }
}

static void test_sigma_inel(FILE *fp,real rho)
{
  int  i;

  fprintf(fp,"Testing the inelastic cross sections table\n");
  for(i=0; (i<500); i++) {
    fprintf(fp,"%3d  %8.3f\n",i,cross_inel(i,rho,NULL));
  }
}

static void test_band_ener(int *seed,FILE *fp)
{
  int i;
  
  for(i=0; (i<500); i++) {
    band_ener(seed,fp,NULL);
  }
}

void init_tables(int nfile,t_filenm fnm[],real rho)
{
  int  seed  = 1993;
  real ekin  = 20;
  real omega = 10;
  
  (void) band_ener(&seed,NULL,opt2fn("-band",nfile,fnm));
  (void) cross_el(ekin,rho,opt2fn("-sigel",nfile,fnm));
  (void) cross_inel(ekin,rho,opt2fn("-sigin",nfile,fnm));
  (void) get_theta_el(ekin,&seed,NULL,opt2fn("-thetael",nfile,fnm));
  (void) get_omega(ekin,&seed,NULL,opt2fn("-eloss",nfile,fnm));
  (void) get_q_inel(ekin,omega,&seed,NULL,opt2fn("-qtrans",nfile,fnm));
}

void test_tables(int *seed,char *fn,real rho)
{
  FILE *fp;
  
  fp = fopen(fn,"w");

  test_omega(fp,seed);
  test_q_inel(fp,seed);
  test_theta_el(fp,seed);
  test_sigma_el(fp,rho);
  test_sigma_inel(fp,rho);
  test_band_ener(seed,fp);
  
  fclose(fp);
}

