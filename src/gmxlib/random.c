/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Grunge ROck MAChoS
 */
static char *SRCID_random_c = "$Id$";

#include <math.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "physics.h"
#include "typedefs.h"
#include "vec.h"
#include "random.h"

real gauss (real am, real sd, int *ig)
     /* we use the boxmuller method  to calculate the gaussian.    */
     /* this method calculates two gaussians at a time, one of     */
     /* which is stored in glgset until the next call to gauss.    */
     /* see: numerical recipes p.716                               */
{
  /* use our own gaussian, or modified boxmuller's ? */
#define STRAND

#ifdef STRAND
  real a;
  int  i;
#else
  static bool gliset=FALSE;
  static real glgset;
  real fac,r,v1,v2;
#endif
  real gval;

#ifdef STRAND
  a = 0;
  for (i=0; (i<12); i++)
    a = a + rando (ig);
  gval = (a-6.0)*sd/6.0+am;
#else
  if (!gliset) {
    do {
      v1 = 2.0*rando(ig)-1.0;
      v2 = 2.0*rando(ig)-1.0;
      r = v1*v1 + v2*v2;
    } while (r >= 1.0);
    fac = sqrt (-2*log(r)/r);
    glgset = v1*fac;
    gval = v2*fac;
  }
  else 
    gval = glgset; /* use the spare gaussian  */
  gliset=!gliset;  /* set the spare indicator */
#endif

  return gval;
}

void low_mspeed(real tempi,int nrdf,int seed,int nat,atom_id a[],
		t_atoms *atoms,rvec v[])
{
  int  i,j,m,ig;
  real boltz,sd;
  real ekin,temp,mass,scal;

  ig=seed;
  boltz=BOLTZ*tempi*nrdf/((real) (DIM*nat));
  ekin=0.0;
  for (i=0; (i<nat); i++) {
    j=a[i];
    mass=atoms->atom[j].m;
    if (mass > 0) {
      sd=sqrt(boltz/mass);
      for (m=0; (m<DIM); m++) {
	v[j][m]=gauss(0.0,sd,&ig);
	ekin+=0.5*mass*v[j][m]*v[j][m];
      }
    }
  }
  temp=(2.0*ekin)/(nrdf*BOLTZ);
  if (temp > 0) {
    scal=tempi/temp;
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

  for(i=0; (i<grp->nr); i++) {
    s=grp->index[i];
    n=grp->index[i+1]-s;
    low_mspeed(tempi[i],nrdf[i],seed,n,&(grp->a[s]),atoms,v);
  }
}

void maxwell_speed(real tempi,int nrdf,int seed,t_atoms *atoms, rvec v[])
{
  atom_id *dummy;
  int     i;

  snew(dummy,atoms->nr);
  for(i=0; (i<atoms->nr); i++)
    dummy[i]=i;
  low_mspeed(tempi,nrdf,seed,atoms->nr,dummy,atoms,v);
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
