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
 * Gnomes, ROck Monsters And Chili Sauce
 */

#include <math.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "physics.h"
#include "typedefs.h"
#include "vec.h"
#include "random.h"

#define GAUSS_NXP 16
/* The size of the data area is 2^GAUSS_NXP */

struct t_gaussdata {
  real  *x;      /* Pointer to the work area */
  int   seed;    /* The random seed */
  int   uselast; /* Is there a saved number we can use? */
  real  last;    /* The possibly saved number */
};


/* Initialize (and warm up) a gaussian random number generator
 * by copying the seed. The routine returns a handle to the
 * new generator.
 */
t_Gaussdata 
init_gauss(int seed)
{
  int size = (1 << GAUSS_NXP);
  int nwarmup = GAUSS_NXP*(size-1);
  int k;
  real tmp;
  t_Gaussdata gaussdata;
  gaussdata=(t_Gaussdata)malloc(sizeof(struct t_gaussdata));
  gaussdata->x=(real *)malloc(size*sizeof(real));
  gaussdata->seed    = seed;
  gaussdata->last    = 0;
  gaussdata->uselast = 0;
  
  for(k=0;k<size;k++)
    gaussdata->x[k]=1;

  for(k=0;k<nwarmup;k++)
    tmp=gauss(gaussdata);

  return gaussdata;
}



/* Return a new gaussian random number with expectation value
 * 0.0 and standard deviation 1.0. This routine is NOT thread-safe
 * for performance reasons - you will either have to do the locking
 * yourself, or better: initialize one generator per thread.
 */
real 
gauss(t_Gaussdata gaussdata)
{
  int i,n1=0,n2;
  int intt;
  int mo;
  int j1;
  int isgn;
  int ne = 31-GAUSS_NXP;

  if(gaussdata->uselast) {
    gaussdata->uselast=0;
    return gaussdata->last;
  } else {
    do {
      for (i=0;i<2;i++) {
	intt=(gaussdata->seed)/127773;
	mo=(gaussdata->seed)-intt*127773;
	j1=2836*intt;
	(gaussdata->seed)=16807*mo-j1;
	if(gaussdata->seed<0)
	  gaussdata->seed+=2147483647;
	if(i==0)
	  n1 = gaussdata->seed >> ne;
      }
      n2 = gaussdata->seed >> ne;
    } while (n1==n2);
    
    isgn=2*(gaussdata->seed & 1)-1;

    gaussdata->x[n1]=isgn*(0.7071067811865475*(gaussdata->x[n1]+gaussdata->x[n2]));
    gaussdata->x[n2]=-gaussdata->x[n1]+isgn*1.414213562373095*gaussdata->x[n2];
    gaussdata->last = gaussdata->x[n2];
    gaussdata->uselast = 1;
    return gaussdata->x[n1];
  }
}



/* Release all the resources used for the generator */
void 
finish_gauss(t_Gaussdata data)
{
  free(data->x);
  free(data);
  data=NULL;
  
  return;
}



void low_mspeed(real tempi,int nrdf,int nat,atom_id a[],
		t_atoms *atoms,rvec v[], t_Gaussdata gaussdata)
{
  int  i,j,m;
  real boltz,sd;
  real ekin,temp,mass,scal;

  boltz=BOLTZ*tempi;
  ekin=0.0;
  for (i=0; (i<nat); i++) {
    j=a[i];
    mass=atoms->atom[j].m;
    if (mass > 0) {
      sd=sqrt(boltz/mass);
      for (m=0; (m<DIM); m++) {
	v[j][m]=sd*gauss(gaussdata);
	ekin+=0.5*mass*v[j][m]*v[j][m];
      }
    }
  }
  temp=(2.0*ekin)/(nrdf*BOLTZ);
  if (temp > 0) {
    scal=sqrt(tempi/temp);
    for(i=0; (i<nat); i++)
      for(m=0; (m<DIM); m++)
	v[a[i]][m]*=scal;
  }
  fprintf(stderr,"Velocities were taken from a Maxwell distribution at %g K\n",
	  tempi);
  if (debug) {
    fprintf(debug,
	    "Velocities were taken from a Maxwell distribution\n"
	    "Initial generated temperature: %12.5e (scaled to: %12.5e)\n",
	    temp,tempi);
  }
}

void grp_maxwell(t_block *grp,real tempi[],int nrdf[],int seed,
		 t_atoms *atoms,rvec v[])
{
  int i,s,n;
  t_Gaussdata gaussdata;
  bool bFirst=TRUE;

  if(bFirst) {
    bFirst=FALSE;
    gaussdata = init_gauss(seed);
  }

  for(i=0; (i<grp->nr); i++) {
    s=grp->index[i];
    n=grp->index[i+1]-s;
    low_mspeed(tempi[i],nrdf[i],n,&(grp->a[s]),atoms,v,gaussdata);
  }
}


void maxwell_speed(real tempi,int nrdf,int seed,t_atoms *atoms, rvec v[])
{
  atom_id *dummy;
  int     i;
  t_Gaussdata gaussdata;
  bool bFirst=TRUE;
  
  if (seed == -1) {
    seed = make_seed();
    fprintf(stderr,"Using random seed %d for generating velocities\n",seed);
  }

  if(bFirst) {
    bFirst=FALSE;
    gaussdata = init_gauss(seed);
  }

  snew(dummy,atoms->nr);
  for(i=0; (i<atoms->nr); i++)
    dummy[i]=i;
  low_mspeed(tempi,nrdf,atoms->nr,dummy,atoms,v,gaussdata);
  sfree(dummy);
}

real calc_cm(FILE *log,int natoms,real mass[],rvec x[],rvec v[],
	     rvec xcm,rvec vcm,rvec acm,matrix L)
{
  rvec dx,a0;
  real tm,m0;
  int  i,m;

  clear_rvec(xcm);
  clear_rvec(vcm);
  clear_rvec(acm);
  tm=0.0;
  for(i=0; (i<natoms); i++) {
    m0=mass[i];
    tm+=m0;
    oprod(x[i],v[i],a0);
    for(m=0; (m<DIM); m++) {
      xcm[m]+=m0*x[i][m]; /* c.o.m. position */
      vcm[m]+=m0*v[i][m]; /* c.o.m. velocity */
      acm[m]+=m0*a0[m];   /* rotational velocity around c.o.m. */
    }
  }
  oprod(xcm,vcm,a0);
  for(m=0; (m<DIM); m++) {
    xcm[m]/=tm;
    vcm[m]/=tm;
    acm[m]-=a0[m]/tm;
  }

#define PVEC(str,v) fprintf(log,\
			    "%s[X]: %10.5e  %s[Y]: %10.5e  %s[Z]: %10.5e\n", \
			   str,v[0],str,v[1],str,v[2])
#ifdef DEBUG
  PVEC("xcm",xcm);
  PVEC("acm",acm);
  PVEC("vcm",vcm);
#endif
  
  clear_mat(L);
  for(i=0; (i<natoms); i++) {
    m0=mass[i];
    for(m=0; (m<DIM); m++)
      dx[m]=x[i][m]-xcm[m];
    L[XX][XX]+=dx[XX]*dx[XX]*m0;
    L[XX][YY]+=dx[XX]*dx[YY]*m0;
    L[XX][ZZ]+=dx[XX]*dx[ZZ]*m0;
    L[YY][YY]+=dx[YY]*dx[YY]*m0;
    L[YY][ZZ]+=dx[YY]*dx[ZZ]*m0;
    L[ZZ][ZZ]+=dx[ZZ]*dx[ZZ]*m0;
  }
#ifdef DEBUG
  PVEC("L-x",L[XX]);
  PVEC("L-y",L[YY]);
  PVEC("L-z",L[ZZ]);
#endif

  return tm;
}

void stop_cm(FILE *log,int natoms,real mass[],rvec x[],rvec v[])
{
  rvec   xcm,vcm,acm;
  tensor L;
  int    i,m;

#ifdef DEBUG
  fprintf(log,"stopping center of mass motion...\n");
#endif
  (void)calc_cm(log,natoms,mass,x,v,xcm,vcm,acm,L);
  
  /* Subtract center of mass velocity */
  for(i=0; (i<natoms); i++)
    for(m=0; (m<DIM); m++)
      v[i][m]-=vcm[m];

#ifdef DEBUG
  (void)calc_cm(log,natoms,mass,x,v,xcm,vcm,acm,L);
#endif
}
