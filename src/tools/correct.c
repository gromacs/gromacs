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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "vec.h"
#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "disco.h"
#include "xtcio.h"
#include "main.h"
#include "random.h"
#include "bondf.h"
#include "names.h"
#include "physics.h"

#define NOT32
#ifdef NOT32
#define myrand(seed) rando(seed)
#else
static unsigned int IA=16807;
static unsigned int IM=2147483647;
static unsigned int IQ=127773;
static unsigned int IR=2836;
static unsigned int MASK=123459876;
static float AM;

void init_rand()
{
  AM = (1.0/IM);
}

float myrand(int *seed)
{
  int   k;
  float ans;

  *seed ^= MASK;
  k      = (*seed)/IQ;
  *seed  = IA*(*seed-k*IQ)-IR*k;
  if (*seed < 0) *seed += IM;
  ans    = AM*(*seed);
  *seed ^= MASK;
  
  return ans;
}
#endif

real mygauss(int *seed,real lb,real ub)
{
  int  i;
  real myg=0.0;
  
  for(i=0; (i<8); i++)
    myg += myrand(seed);
    
  return lb+0.125*(ub-lb);
}

void randomize_list(int n,int ip[],int *seed)
{
  int i,iran,temp;
  
  for(i=0; (i<n); i++) {
    iran     = ((int) (myrand(seed)*n)) % n;
    temp     = ip[i];
    ip[i]    = ip[iran];
    ip[iran] = temp;
  }
}

void reset_list(int n,int ip[])
{
  int i;
  
  for(i=0; (i<n); i++) 
    ip[i]    = i;
}

int update_list(t_dist d[],int tag[],int natom,bool bViol[],int ip[])
{
  int i,j,j0,j1,nchk;
  
  nchk = 0;
  for(i=0; (i<natom); i++) {
    j0 = tag[i];
    j1 = tag[i+1];
    if (bViol[i]) {
      for(j=j0; (j<j1); j++) {
	if (d[j].ai != i)
	  gmx_fatal(FARGS,"Tags or distances inconsistent: "
		      "i=%d, j=%d, ai=%d, aj=%d, j0=%d, j1=%d",
		      i,j,d[j].ai,d[j].aj,j0,j1);
	ip[nchk] = j;
	nchk++;
      }
    }
    else {
      for(j=j0; (j<j1); j++) {
	if (d[j].ai != i)
	  gmx_fatal(FARGS,"Tags or distances inconsistent: "
		      "i=%d, j=%d, ai=%d, aj=%d, j0=%d, j1=%d",
		      i,j,d[j].ai,d[j].aj,j0,j1);
	if (bViol[d[j].aj]) {
	  ip[nchk] = j;
	  nchk++;
	}
      }
    }
  }
  return nchk;
}

bool mirror_xx(int n1,int n2,int n3,int n4,rvec x[],bool bViol[],real xi)
{
  /* Check the impropers, and fix if necessary */
  /* adapt the coordinates of atom n4 such that it is mirrorred
   * with respect to the plane defined by n1,n2,n3 (n1 origin)
   */
  rvec v1,v2,x4old,n,d14;
  real lambda,nsq,dn;
  int  m;
  
  copy_rvec(x[n4],x4old);
  
  /* v1 and v2 span the plane */
  rvec_sub(x[n2],x[n1],v1);
  rvec_sub(x[n3],x[n1],v2);
  rvec_sub(x[n1],x[n4],d14);

  /* Take the normal vector to the plane */
  oprod(v1,v2,n);
  nsq = iprod(n,n);
  
  /* Check for x[n1], x[n2], x[n3] on a line, if so come back next iteration */
  if (nsq == 0.0)
    return FALSE;
    
  /* Find the point in which the vector parallel to n thru x[n4]
   * intersects with the plane:
   * by solving (x[n4]+lambda*n) . n = x[n1] . n 
   * Compute lambda 
   */
  lambda = iprod(d14,n)/nsq;

  /* Check on the improper: If it is really far off than we should not
   * try to correct if, as this may only make it worse. The limits below
   * indicate a range that when mirrored will become a normal improper
   * again.
   */
  if ((xi < -10) && (xi > -60)) {
    /* The fourth atom is at the wrong side of the plane, fix it. */
  
    /* Mirror all four atoms around the plane parallel to 123, going thru
     * the center of mass of the four atoms. It is easy to see that there
     * is no center-of-mass displacement this way.
     */
    for(m=0; (m<DIM); m++) {
      dn = 0.5*lambda*n[m];
      x[n1][m] -= dn;
      x[n2][m] -= dn;
      x[n3][m] -= dn;
      x[n4][m] += 3*dn;
    }
    if (debug) {
      fprintf(debug,"improper: %d %d %d %d, xi = %g deg.\n",
	      n1+1,n2+1,n3+1,n4+1,xi);
      fprintf(debug,"lambda = %g, n = (%g,%g,%g)\n",lambda,n[XX],n[YY],n[ZZ]);
      fprintf(debug,"x4old = (%g,%g,%g) new = (%g,%g,%g)\n",
	      x4old[XX],x4old[YY],x4old[ZZ],x[n4][XX],x[n4][YY],x[n4][ZZ]);
      bViol[n1] = TRUE;
      bViol[n2] = TRUE;
      bViol[n3] = TRUE;
      bViol[n4] = TRUE;
    }
    return TRUE;
  }
  return FALSE;
}

void compute_dihs(FILE *log,int nq,t_quadruple q[],rvec x[],real phi[],
		  matrix box)
{
  int    i,t1,t2,t3;
  rvec   r_ij,r_kj,r_kl,m,n;
  real   cos_phi,sign;
  
  for(i=0; (i<nq); i++)
    phi[i] = RAD2DEG*dih_angle(x[q[i].ai],x[q[i].aj],x[q[i].ak],x[q[i].al],
			       NULL,
			       r_ij,r_kj,r_kl,m,n,&cos_phi,&sign,
			       &t1,&t2,&t3);
}

int check_impropers(FILE *log,t_correct *c,int natom,rvec x[],matrix box)
		    
{
  int  i,nmirror,nwrong;
  bool bFlip = FALSE;
  
  nmirror = 0;
  if (c->bChiral) {
    compute_dihs(log,c->nimp,c->imp,x,c->idih,box);
    nwrong = 0;
    for(i=0; (i<c->nimp); i++)
      if (c->idih[i] < 0)
	nwrong++;
    if (nwrong > c->nimp/2) {
      if (debug)
	fprintf(debug,"Flipping entire molecule\n");
      for(i=0; (i<natom); i++)
	svmul(-1.0,x[i],x[i]);
    }
    for(i=0; (i<c->nimp); i++)
      if (mirror_xx(c->imp[i].ai,c->imp[i].aj,c->imp[i].ak,c->imp[i].al,
		    x,c->bViol,c->idih[i]))
	nmirror++;
    if (debug && (nmirror > 0))
      fprintf(debug,"mirrored %d improper dihedrals\n",nmirror);
  }
  if (bFlip)
    nmirror = c->nimp - nmirror;
    
  return nmirror;
}

int check_cis_trans(FILE *log,t_correct *c,rvec x[],matrix box)
{
  int  i,nswap;
  rvec temp;
  
  nswap = 0;
  if (c->bPep) {
    compute_dihs(log,c->npep,c->pepbond,x,c->omega,box);
    for(i=0; (i<c->npep); i++) {
      /* If one of the peptide bonds is cis, then swap H and Ca */
      if ((c->omega[i] < -90) || (c->omega[i] > 90)) {
	copy_rvec(x[c->pepbond[i].al],temp);
	copy_rvec(x[c->pepbond[i].am],x[c->pepbond[i].al]);
	copy_rvec(temp,x[c->pepbond[i].am]);
	/* These are now violated! */
	c->bViol[c->pepbond[i].al] = TRUE;
	c->bViol[c->pepbond[i].am] = TRUE;
	nswap++;
      }
    }
    if (debug && (nswap > 0))
      fprintf(debug,"swapped %d peptide bonds\n",nswap);
  }
  return nswap;
}

bool do_1steep(int  ai,      int aj,       rvec dx,
	       real len,     real lb,      real ub,     real wi,
	       real *dev,    real *maxdev, real *ener, rvec fptr[],
	       bool bViol[])
{
  int  m;
  real ndev,wj,guess,ddx,delta;
  
  /* Violated bounds ? */
  if ((len < lb) || (len > ub)) {
    /* Update counters */
    if (len < lb) {
      ndev  = lb-len;
      guess = lb;
    }
    else {
      ndev  = len-ub;
      guess = ub;
    }
    
    *dev    += ndev;
    *ener   += (ndev*ndev);
    *maxdev  = max(*maxdev,ndev);
    
    /* Set violation flag */
    bViol[ai] = TRUE;
    bViol[aj] = TRUE;
    
    /* Correct non-bondeds only every N steps */
    wj  = 1.0-wi;
    
    /* Calculate the force, that is the deviation from the nearest
     * bound. This is the steepest descent minimizaton part.
     */
    delta = (len-guess)/len;
    for(m=0; (m<DIM); m++) {
      ddx = delta*dx[m];
      fptr[ai][m] += wi*ddx;
      fptr[aj][m] -= wj*ddx;
    }
    return TRUE;
  }
  else
    return FALSE;
}

bool do_1shake(int  cons_type,
	       int  ai,      int aj,       rvec dx,     
	       real len,     real lb,      real ub,     real wi,
	       real *dev,    real *maxdev, real *ener,  rvec xptr[],
	       bool bViol[], int  *seed,   bool bLowerOnly)
{
  int  m;
  real ndev,delta,ddx,wj,guess;
  
  /* Violated bounds ? */
  if ((len < lb) || ((cons_type != edcNONBOND) && (len > ub))) {
    /* Update counters */
    if (cons_type == edcNONBOND) 
      ndev = lb-len;
    else {
      if (len < lb) 
	ndev = lb-len;
      else
	ndev = len-ub;
    }
    *dev    += ndev;
    *ener   += (ndev*ndev);
    *maxdev  = max(*maxdev,ndev);
    
    /* Set violation flag */
    bViol[ai] = TRUE;
    bViol[aj] = TRUE;
	
    wj  = 1.0-wi;
    if (cons_type == edcBOND)
      guess = mygauss(seed,lb,ub);
    else
      guess = lb + myrand(seed)*(ub-lb);
    
    delta = (len-guess)/len;
    for(m=0; (m<DIM); m++) {
      ddx = delta*dx[m];
      xptr[ai][m] += wi*ddx;
      xptr[aj][m] -= wj*ddx;
    }
    return TRUE;
  }
  else
    return FALSE;
}

int shake_coords(FILE *log,bool bVerbose,
		 int nstruct,int natom,rvec xref[],
		 rvec x[],int *seed,matrix box,t_correct *c,int *niter)
{
  static bool bFirst=TRUE;
  static rvec *force[2],*xtry[2],*xwr;
  int    nit,nchk,low,oldlow;
  int    i,j,k,l,m,ai,aj,status,nprint;
  int    nviol[edcNR],nvtot=0,nswap,nmirror,cons_type;
  bool   bShort,bConverged,bNB,bPrint,bLowDev,bExplicit,bMyNB,bLowerOnly;
  bool   bVNow;
  bool   *bViol;
  real   lb,ub,lbfac,wi,len,dev,maxdev,len2;
  real   step,maxforce=0,rstep,step0;
  real   ener[2];
  int    cur=0;
#define next 1-cur
  t_dist *d;
  int    *ip;
  rvec   dx,xcm,*xptr,*fptr;
  char   buf[32];
  
  /* Initiate local variables */
  bViol      = c->bViol;
  d          = c->d;
  ip         = c->ip;
  nchk       = c->ndist;
  low        = c->ndist;
  oldlow     = c->ndist;
  bExplicit  = c->bExplicit;
  bConverged = FALSE;
  bLowDev    = FALSE;
  bLowerOnly = c->bLowerOnly;
  step0      = 0.1;
  step       = step0;
  nprint     = 0;
  ener[cur]  = 1e30;
  status     = 0;
  
  debug_gmx();  
  if (bFirst) {
    snew(force[0],natom);
    snew(force[1],natom);
    snew(xtry[0],natom);
    snew(xtry[1],natom);
    snew(xwr,natom);
#ifndef NOT32
    init_rand();
#endif
    bFirst = FALSE;
  }
  if (c->bDump) {  
    sprintf(buf,"test%d.xtc",nstruct);
    status     = open_xtc(buf,"w");
    center_in_box(natom,x,box,xwr);
    write_xtc(status,natom,0,0,box,xwr,1000.0);
  }
  for(i=0; (i<natom); i++) {
    copy_rvec(x[i],xtry[cur][i]);
    copy_rvec(x[i],xtry[next][i]);
  }
  xptr = xtry[next];
  debug_gmx();
  
  if (c->bRanlistFirst)
    randomize_list(nchk,ip,seed);
  
  for(nit=0; ((nit < c->maxnit) && !bConverged); nit++) {
    /* Reset violation counters */
    for(i=0; (i<natom); i++)
      bViol[i] = FALSE;

    /* Set progress indicators */
    nviol[edcNONBOND] = 0;
    nviol[edcBOND]    = 0;
    nviol[edcDISRE]   = 0;
    dev        = 0.0;
    maxdev     = 0.0;
    ener[next] = 0.0;

    /* Nonbonded growing factor */
    lbfac      = (c->ngrow == 0) ? 1.0 : (min(1.0,(real)nit/(real)c->ngrow));
        
    /* Set step dependent flags */
    bNB        = (c->nbcheck != 0) && (nit % c->nbcheck)  == 0;
    bPrint     = (c->nstprint!= 0) && (nit % c->nstprint) == 0;
    bExplicit  = bLowDev && c->bExplicit;
    
    if (c->nstranlist && ((nit % c->nstranlist) == 0))
      randomize_list(nchk,ip,seed);
      
    xptr = xtry[next];
    if (bExplicit) {
      fptr = force[next];
      for(i=0; (i<natom); i++)
	clear_rvec(fptr[i]);
    }
    else 
      fptr = NULL;
    
    /* Loop over distance list */
    for(l=0; (l<nchk); l++) {
      /* Unpack lists */
      k         = ip[l];
      ai        = d[k].ai;
      aj        = d[k].aj;
      lb        = d[k].lb;
      ub        = d[k].ub;
      wi        = d[k].wi;
      cons_type = d[k].cons_type;
      
      if (cons_type == edcNONBOND) {
      	lb *= lbfac;
	bMyNB = bNB;
      }
      else
	bMyNB = FALSE;
      
      if (bMyNB || (cons_type != edcNONBOND)) {
	/* Compute distance */
	len2 = 0.0;
	for(m=0; (m<DIM); m++) {
	  dx[m] = xptr[aj][m] - xptr[ai][m];
	  len2 += dx[m]*dx[m];
	}
	len = sqrt(len2);
	
	if (bExplicit)
	  bVNow = do_1steep(ai,aj,dx,len,lb,ub,wi,&dev,&maxdev,&ener[next],
			    fptr,bViol);
	else {
	  bVNow = do_1shake(cons_type,ai,aj,
			    dx,len,lb,ub,wi,&dev,&maxdev,&ener[next],
			    xptr,bViol,seed,bLowerOnly);
#ifdef HEAVY
	  if (bVNow && c->bDump && debug) {
	    center_in_box(natom,xptr,box,xwr);
	    write_xtc(status,natom,nit,(real)l,box,xwr,1000.0);
	  }
#endif
	}
	if (bVNow)
	  nviol[cons_type]++;
      }
    }
    if (bExplicit) {
      if (ener[next] < ener[cur]) {
	cur = next;
	maxforce = 0;
	for(i=0; (i<natom); i++)
	  for(m=0; (m<DIM); m++)
	    maxforce = max(maxforce,fabs(fptr[i][m]));
	/* Increase step again */
	step  = min(step0,2*step);
	rstep = step/maxforce;
	for(i=0; (i<natom); i++)
	  for(m=0; (m<DIM); m++)
	    xtry[next][i][m] = xtry[cur][i][m] + fptr[i][m]*rstep;
	xptr = xtry[next];
      }
      else {
	/* Decrease step size and try again */
	step  = step/2;
	rstep = step/maxforce;
	for(i=0; (i<natom); i++)
	  for(m=0; (m<DIM); m++)
	    xtry[next][i][m] = xtry[cur][i][m] + fptr[i][m]*rstep;
      }
    }
    if (c->bDump && bPrint) {
      center_in_box(natom,xptr,box,xwr);
      write_xtc(status,natom,nit,(real)nit+1,box,xwr,1000.0);
    }
    /* Check Impropers */
    nmirror = check_impropers(bVerbose ? log : NULL,c,natom,xptr,box);
    
    /* Check peptide bonds */
    nswap   = check_cis_trans(bVerbose ? log : NULL,c,xptr,box);
    
    nvtot   = (nviol[edcNONBOND] + nviol[edcBOND] + nviol[edcDISRE] + 
	       nmirror + nswap);
    if (bPrint) {
      if ((nprint % 20) == 0) 
	fprintf(stderr,"\n%5s%7s%7s%7s%7s%7s%7s%7s %10s %10s\n",
		"nit","nchk","nimp","npep",
		edc_names[edcNONBOND],edc_names[edcBOND],edc_names[edcDISRE],
		"nvtot","maxdev",
		bExplicit ? "ener A^2" : "<dev>");
      fprintf(stderr,"%5d%7d%7d%7d%7d%7d%7d%7d %10.4e %10.4e\n",
	      nit,nchk,nmirror,nswap,
	      nviol[edcNONBOND],nviol[edcBOND],nviol[edcDISRE],nvtot,maxdev,
	      bExplicit ? ener[cur] : dev/natom);
      nprint++;
    }
    bLowDev    = (bLowDev || (dev/natom < c->lodev));
    
    /* Update neighbour list */
#define UPDL
#ifdef UPDL
    nchk = update_list(c->d,c->tag,natom,bViol,ip);
    
    bConverged = (nchk == 0);
#else
    reset_list(nchk,ip);
    
    bConverged = (nvtot == 0);
#endif
  }
  if (c->bDump)
    close_xtc(status);

  /* Copy final structure back to input array */
  for(i=0; (i<natom); i++) 
    copy_rvec(xptr[i],x[i]);
  
  fprintf(log,"struct: %5d, box: %8.3f  %8.3f  %8.3f\n",
	  nstruct,box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
  fprintf(log,"Converged: %3s nit: %5d  <dev>: %8.4f  maxdev: %8.4f\n",
	  yesno_names[bConverged],nit,dev/natom,maxdev);

  *niter = nit;
	    
  return nvtot;
}
    
int quick_check(FILE *log,int natom,rvec x[],matrix box,t_correct *c)
{
  int  i,j,k,m,ai,aj,nviol;
  real lb,ub,len,len2,dx;

  if (log) 
    fprintf(log,"Double checking final structure\n");
  nviol = 0;
  for(k=0; (k<c->ndist); k++) {
    /* Unpack lists */
    ai    = c->d[k].ai;
    aj    = c->d[k].aj;
    lb    = c->d[k].lb;
    ub    = c->d[k].ub;

    /* Compute distance */
    len2 = 0.0;
    for(m=0; (m<DIM); m++) {
      dx    = x[aj][m] - x[ai][m];
      len2 += dx*dx;
    }
    len = sqrt(len2);
    if ((len < lb) || ((ub > 0) && (len > ub))) {
      if (debug)
	fprintf(debug,"VIOL: ai=%4d aj=%4d lb=%8.4f len=%8.4f ub=%8.4f\n",
		ai,aj,lb,len,ub);
      nviol++;
    }
  }
  /* Check Impropers */
  /* nmirror = check_impropers(log,c,natom,x,box);*/
  
  /* Check peptide bonds */
  /* nswap   = check_cis_trans(log,c,x,box);*/
  
  return nviol;
}
