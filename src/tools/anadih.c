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
 * Gnomes, ROck Monsters And Chili Sauce
 */
static char *SRCID_anadih_c = "$Id$";

#include <math.h>
#include <stdio.h>
#include "physics.h"
#include "smalloc.h"
#include "macros.h"
#include "txtdump.h"
#include "bondf.h"
#include "xvgr.h"
#include "typedefs.h"
#include "gstat.h"
#include "lutab.h"

static int calc_RBbin(real phi)
{
  static const real r30  = M_PI/6.0;
  static const real r90  = M_PI/2.0;
  static const real r150 = M_PI*5.0/6.0;
  
  if ((phi < r30) && (phi > -r30))
    return 1;
  else if ((phi > -r150) && (phi < -r90))
    return 2;
  else if ((phi < r150) && (phi > r90))
    return 3;
  return 0;
}

static int calc_Nbin(real phi)
{
  static const real r30  =  30*DEG2RAD;
  static const real r90  =  90*DEG2RAD;
  static const real r150 = 150*DEG2RAD;
  static const real r210 = 210*DEG2RAD;
  static const real r270 = 270*DEG2RAD;
  static const real r330 = 330*DEG2RAD;
  static const real r360 = 360*DEG2RAD;
  
  if (phi < 0)
    phi += r360;
    
  if ((phi > r30) && (phi < r90))
    return 1;
  else if ((phi > r150) && (phi < r210))
    return 2;
  else if ((phi > r270) && (phi < r330))
    return 3;
  return 0;
}

void ana_dih_trans(char *fn_trans,char *fn_histo,
		   real **dih,int nframes,int nangles,
		   char *grpname,real t0,real dt,bool bRb)
{
  FILE *fp;
  int  *tr_f,*tr_h;
  char title[256];
  int  i,j,ntrans;
  int  cur_bin,new_bin;
  real ttime,tt,mind, maxd, prev;
  int  (*calc_bin)(real);
  
  /* Analysis of dihedral transitions */
  fprintf(stderr,"Now calculating transitions...\n");

  if (bRb)
    calc_bin=calc_RBbin;
  else
    calc_bin=calc_Nbin;
    
  snew(tr_h,nangles);
  snew(tr_f,nframes);
    
  /* dih[i][j] is the dihedral angle i in frame j  */
  ntrans = 0;
  for (i=0; (i<nangles); i++)
  {
    /*#define OLDIE*/
#ifdef OLDIE
    mind = maxd = prev = dih[i][0]; 
#else
    cur_bin = calc_bin(dih[i][0]);
#endif    
    for (j=1; (j<nframes); j++)
    {
      new_bin = calc_bin(dih[i][j]);
#ifndef OLDIE
      if (cur_bin == 0)
	cur_bin=new_bin;
      else if ((new_bin != 0) && (cur_bin != new_bin)) {
	cur_bin = new_bin;
	tr_f[j]++;
	tr_h[i]++;
	ntrans++;
      }
#else
      /* why is all this md rubbish periodic? Remove 360 degree periodicity */
      if ( (dih[i][j] - prev) > M_PI)
	dih[i][j] -= 2*M_PI;
      else if ( (dih[i][j] - prev) < -M_PI)
	dih[i][j] += 2*M_PI;

      prev = dih[i][j];
       
      mind = min(mind, dih[i][j]);
      maxd = max(maxd, dih[i][j]);
      if ( (maxd - mind) > 2*M_PI/3)    /* or 120 degrees, assuming       */
      {                                 /* multiplicity 3. Not so general.*/	
	tr_f[j]++;
	tr_h[i]++;
	maxd = mind = dih[i][j];        /* get ready for next transition  */
	ntrans++;
      }
#endif
    }
  }
  fprintf(stderr,"Total number of transitions: %10d\n",ntrans);
  if (ntrans > 0) {
    ttime = (dt*nframes*nangles)/ntrans;
    fprintf(stderr,"Time between transitions:    %10.3f ps\n",ttime);
  }
    
  sprintf(title,"Number of transitions: %s",grpname);
  fp=xvgropen(fn_trans,title,"Time (ps)","# transitions/timeframe");
  for(j=0; (j<nframes); j++) {
    tt = t0+j*dt;
    fprintf(fp,"%10.3f  %10d\n",tt,tr_f[j]);
  }
  ffclose(fp);

  /* Compute histogram from # transitions per dihedral */
  /* Use old array */
  for(j=0; (j<nframes); j++)
    tr_f[j]=0;
  for(i=0; (i<nangles); i++)
    tr_f[tr_h[i]]++;
  for(j=nframes; ((tr_f[j-1] == 0) && (j>0)); j--)
    ;
  
  ttime = dt*nframes;
  sprintf(title,"Transition time: %s",grpname);
  fp=xvgropen(fn_histo,title,"Time (ps)","#");
  for(i=j-1; (i>0); i--) {
    if (tr_f[i] != 0)
      fprintf(fp,"%10.3f  %10d\n",ttime/i,tr_f[i]);
  }
  ffclose(fp);
  
  sfree(tr_f);
  sfree(tr_h);
}
    
void calc_distribution_props(int nh,int histo[],real start,
			     int nkkk, t_karplus kkk[],
			     real *S2)
{
  real d,dc,ds,c1,c2,tdc,tds;
  real fac,ang,invth;
  int  i,j,th;
  
  if (nh == 0)
    fatal_error(0,"No points in histogram (%s, %d)",__FILE__,__LINE__);
  fac = 2*M_PI/nh;
  
  /* Compute normalisation factor */
  th=0;
  for(j=0; (j<nh); j++) 
    th+=histo[j];
  invth=1.0/th;
  
  for(i=0; (i<nkkk); i++) 
    kkk[i].Jc = 0;
  
  tdc=0,tds=0;
  for(j=0; (j<nh); j++) {
    d    = invth*histo[j];
    ang  = j*fac-start;
    c1   = cos(ang);
    c2   = c1*c1;
    dc   = d*c1;
    ds   = d*sin(ang);
    tdc += dc;
    tds += ds;
    for(i=0; (i<nkkk); i++) {
      c1   = cos(ang+kkk[i].offset);
      c2   = c1*c1;
      kkk[i].Jc += d*(kkk[i].A*c2 + kkk[i].B*c1 + kkk[i].C);
    }
  }
  *S2 = tdc*tdc+tds*tds;
}

static void calc_angles(FILE *log,matrix box,
			int n3,atom_id index[],real ang[],rvec x_s[])
{
  int  i,ix;
  rvec r_ij,r_kj;
  real costh;
  
  for(i=ix=0; (ix<n3); i++,ix+=3) 
    ang[i]=bond_angle(log,box,x_s[index[ix]],x_s[index[ix+1]],x_s[index[ix+2]],
		      r_ij,r_kj,&costh);
  if (debug) {
    fprintf(debug,"Angle[0]=%g, costh=%g, index0 = %u, %u, %u\n",
	    ang[0],costh,index[0],index[1],index[2]);
    pr_rvec(debug,0,"rij",r_ij,DIM);
    pr_rvec(debug,0,"rkj",r_kj,DIM);
    pr_rvecs(debug,0,"box",box,DIM);
  }
}

static real calc_fraction(real angles[], int nangles)
{
  int i;
  real trans = 0, gauche = 0;
  real angle;

  for (i = 0; i < nangles; i++)
  {
    angle = angles[i] * RAD2DEG;

    if (angle > -45 && angle < 45)
      trans += 1.0;
    else if (angle > 90 && angle < 150)
      gauche += 1.0;
    else if (angle < -90 && angle > -150)
      gauche += 1.0;
  }
  if (trans+gauche > 0)
    return trans/(trans+gauche);
  else 
    return 0;
}

static void calc_dihs(FILE *log,matrix box,
		      int n4,atom_id index[],real ang[],rvec x_s[])
{
  int  i,ix;
  rvec r_ij,r_kj,r_kl,m,n;
  real cos_phi,sign,aaa;
  
  for(i=ix=0; (ix<n4); i++,ix+=4) {
    aaa=dih_angle(log,box,
		  x_s[index[ix]],x_s[index[ix+1]],x_s[index[ix+2]],
		  x_s[index[ix+3]],
		  r_ij,r_kj,r_kl,m,n,
		  &cos_phi,&sign);
    ang[i]=aaa;  /* not taking into account ryckaert bellemans yet */
  }
}

void make_histo(FILE *log,
		int ndata,real data[],int npoints,int histo[],
		real minx,real maxx)
{
  double dx;
  int    i,ind;
  
  if (minx == maxx) {
    minx=maxx=data[0];
    for(i=1; (i<ndata); i++) {
      minx=min(minx,data[i]);
      maxx=max(maxx,data[i]);
    }
    fprintf(log,"Min data: %10g  Max data: %10g\n",minx,maxx);
  }
  dx=(double)npoints/(maxx-minx);
  if (debug)
    fprintf(debug,
	    "Histogramming: ndata=%d, nhisto=%d, minx=%g,maxx=%g,dx=%g\n",
	    ndata,npoints,minx,maxx,dx);
  for(i=0; (i<ndata); i++) {
    ind=(data[i]-minx)*dx;
    if ((ind >= 0) && (ind < npoints))
      histo[ind]++;
    else
      fprintf(log,"index = %d, data[%d] = %g\n",ind,i,data[i]);
  }
}

void normalize_histo(int npoints,int histo[],real dx,real normhisto[])
{
  int    i;
  double d,fac;
  
  d=0;
  for(i=0; (i<npoints); i++)
    d+=dx*histo[i];
  if (d==0) {
    fprintf(stderr,"Empty histogram!\n");
    return;
  }
  fac=1.0/d;
  for(i=0; (i<npoints); i++)
    normhisto[i]=fac*histo[i];
}

void read_ang_dih(char *trj_fn,char *tpb_fn,
		  bool bAngles,bool bSaveAll,bool bRb,
		  int maxangstat,int angstat[],
		  int *nframes,real **time,
		  int isize,atom_id index[],
		  real **trans_frac,
		  real **aver_angle,
		  real *dih[])
{
  t_topology *top;
  int        i,angind,status,natoms,nat,total,teller,nangles,nat_trj,n_alloc;
  real       t,fraction,pifac,aa,angle;
  real       *angles[2];
  matrix     box;
  rvec       *x,*x_s;
  int        cur=0;
#define prev (1-cur)

  /* Initiate lookup table for sqrt calculations */
  init_lookup_table(stdout);

  /* Read topology */    
  top     = read_top(tpb_fn);
  natoms  = top->atoms.nr;
  nat_trj = read_first_x(&status,trj_fn,&t,&x,box);
  
  /* Check for consistency of topology and trajectory */
  if (natoms < nat_trj)
    fprintf(stderr,"WARNING! Topology has fewer atoms than trajectory\n");
    
  snew(x_s,nat_trj);
  
  if (bAngles) {
    nangles=isize/3;
    pifac=M_PI;
  }
  else {
    nangles=isize/4;
    pifac=2.0*M_PI;
  }
  snew(angles[cur],nangles);
  snew(angles[prev],nangles);
  
  /* Start the loop over frames */
  total       = 0;
  teller      = 0;
  n_alloc     = 0;
  *time       = NULL;
  *trans_frac = NULL;
  *aver_angle = NULL;

  do {
    if (teller >= n_alloc) {
      n_alloc+=100;
      if (bSaveAll)
	for (i=0; (i<nangles); i++)
	  srenew(dih[i],n_alloc);
      srenew(*time,n_alloc);
      srenew(*trans_frac,n_alloc);
      srenew(*aver_angle,n_alloc);
    }

    (*time)[teller] = t;

    rm_pbc(&(top->idef),nat_trj,box,x,x_s);
    
    if (bAngles)
      calc_angles(stdout,box,isize,index,angles[cur],x_s);
    else {
      calc_dihs(stdout,box,isize,index,angles[cur],x_s);

      /* Trans fraction */
      fraction = calc_fraction(angles[cur], nangles);
      (*trans_frac)[teller] = fraction;
      
      /* Change Ryckaert-Bellemans dihedrals to polymer convention */
      if (bRb)
	for(i=0; (i<nangles); i++)
	  if (angles[cur][i] <= 0.0) 
	    angles[cur][i] += M_PI;
	  else
	    angles[cur][i] -= M_PI;
      
      /* Periodicity in dihedral space... */
      if (teller > 1) {
	for(i=0; (i<nangles); i++) {
	  while (angles[cur][i] <= angles[prev][i] - M_PI)
	    angles[cur][i]+=2*M_PI;
	  while (angles[cur][i] > angles[prev][i] + M_PI)
	    angles[cur][i]-=2*M_PI;
	}
      }
    }

    /* Average angles */      
    aa=0;
    for(i=0; (i<nangles); i++) {
      aa=aa+angles[cur][i];
      
      /* angle in rad / 2Pi * max determines bin. bins go from 0 to maxangstat,
	 even though scale goes from -pi to pi (dihedral) or -pi/2 to pi/2 
	 (angle) Basically: translate the x-axis by Pi. Translate it back by 
	 -Pi when plotting.
       */

      angle = angles[cur][i];
      if (!bAngles) {
	while (angle < -M_PI) 
	  angle += 2*M_PI;
	while (angle >= M_PI) 
	  angle -= 2*M_PI;
	  
	angle+=M_PI;
      }
      
      /* Update the distribution histogram */
      angind = (int) ((angle*maxangstat)/pifac + 0.5);
      if (angind==maxangstat)
	angind=0;
      if ( (angind < 0) || (angind >= maxangstat) )
	/* this will never happen */
	fatal_error(0,"angle (%f) index out of range (0..%d) : %d\n",
		    angle,maxangstat,angind);
      
      angstat[angind]++;
      if (angind==maxangstat)
	fprintf(stderr,"angle %d fr %d = %g\n",i,cur,angle);
      
      total++;
    }
    
    /* average over all angles */
    (*aver_angle)[teller] = (aa/nangles);  
    
    /* this copies all current dih. angles to dih[i], teller is frame */
    if (bSaveAll) 
      for (i = 0; i < nangles; i++)
	dih[i][teller] = angles[cur][i];
      
    /* Swap buffers */
    cur=prev;
    
    /* Increment loop counter */
    teller++;
  } while (read_next_x(status,&t,nat_trj,x,box));  
  close_trj(status); 
  
  sfree(x);
  sfree(x_s);
  sfree(angles[cur]);
  sfree(angles[prev]);
  
  *nframes=teller;
}

