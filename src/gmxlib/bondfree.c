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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_bondfree_c = "$Id$";

#include <math.h>
#include "assert.h"
#include "physics.h"
#include "vec.h"
#include "vveclib.h"
#include "maths.h"
#include "txtdump.h"
#include "bondf.h"
#include "smalloc.h"
#include "pbc.h"
#include "ns.h"
#include "macros.h"
#include "names.h"
#include "fatal.h"
#include "mshift.h"
#include "main.h"
#include "disre.h"

static bool bPBC=FALSE;

void pbc_rvec_sub(rvec xi,rvec xj,rvec dx)
{
  if (bPBC)
    pbc_dx(xi,xj,dx);
  else
    rvec_sub(xi,xj,dx);
}

void calc_bonds(FILE *log,t_commrec *cr,t_idef *idef,
		rvec x_s[],rvec f[],
		t_forcerec *fr,t_graph *g,
		real epot[],t_nrnb *nrnb,
		matrix box,real lambda,
		t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  static bool bFirst=TRUE;
  int    ftype,nbonds,ind,nat;

  if (bFirst) {
    bPBC   = (getenv("GMXFULLPBC") != NULL);
    if (bPBC)
      fprintf(log,"Full PBC calculation = %s\n",bool_names[bPBC]);
    bFirst = FALSE;
#ifdef DEBUG
    p_graph(debug,"Bondage is fun",g);
#endif
  }
  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (interaction_function[ftype].flags & IF_BOND) {
      nbonds=idef->il[ftype].nr;
      if (nbonds > 0) {
	epot[ftype]+=
	  interaction_function[ftype].ifunc(nbonds,idef->il[ftype].iatoms,
					    idef->iparams,x_s,f,fr,g,box,
					    lambda,&epot[F_DVDL],
					    md,ngrp,egnb,egcoul);
	ind = interaction_function[ftype].nrnb_ind;
	nat = interaction_function[ftype].nratoms+1;
	if (ind != -1)
	  inc_nrnb(nrnb,ind,nbonds/nat);
      }
    }
  }
}

/*
 * Morse potential bond by Frank Everdij
 *
 * Three parameters needed:
 *
 * b0 = equilibrium distance in nm
 * be = beta in nm^-1 (actually, it's nu_e*Sqrt(2*pi*pi*mu/D_e))
 * cb = well depth in kJ/mol
 *
 * Note: the potential is referenced to be +cb at infinite separation
 *       and zero at the equilibrium distance!
 */

real morsebonds(int nbonds,
                t_iatom forceatoms[],t_iparams forceparams[],
                rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
                matrix box,real lambda,real *dvdl,
                t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  const real one=1.0;
  const real two=2.0;
  real  dr,dr2,temp,omtemp,cbomtemp,fbond,vbond,fij,b0,be,cb,vtot;
  rvec  dx;
  int   i,m,ki,type,ai,aj;
  ivec  dt;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    
    b0   = forceparams[type].morse.b0;
    be   = forceparams[type].morse.beta;
    cb   = forceparams[type].morse.cb;

    pbc_rvec_sub(x[ai],x[aj],dx);                   /*   3          */
    dr2  = iprod(dx,dx);                            /*   5          */
    dr   = sqrt(dr2);                               /*  10          */
    temp = exp(-be*(dr-b0));                        /*  12          */
    
    if (temp == one)
      continue;

    omtemp   = one-temp;                               /*   1          */
    cbomtemp = cb*omtemp;                              /*   1          */
    vbond    = cbomtemp*omtemp;                        /*   1          */
    fbond    = -two*be*temp*cbomtemp*invsqrt(dr2);      /*   9          */
    vtot    += vbond;       /* 1 */
    
    ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
    ki=IVEC2IS(dt);

    for (m=0; (m<DIM); m++) {                          /*  15          */
      fij=fbond*dx[m];
      f[ai][m]+=fij;
      f[aj][m]-=fij;
      fr->fshift[ki][m]+=fij;
      fr->fshift[CENTRAL][m]-=fij;
    }
  }                                           /*  58 TOTAL    */
  return vtot;
}

real cubicbonds(int nbonds,
                t_iatom forceatoms[],t_iparams forceparams[],
                rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
                matrix box,real lambda,real *dvdl,
                t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  const real three = 3.0;
  const real two   = 2.0;
  real  kb,b0,kcub;
  real  dr,dr2,dist,kdist,kdist2,fbond,vbond,fij,vtot;
  rvec  dx;
  int   i,m,ki,type,ai,aj;
  ivec  dt;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    
    b0   = forceparams[type].cubic.b0;
    kb   = forceparams[type].cubic.kb;
    kcub = forceparams[type].cubic.kcub;

    pbc_rvec_sub(x[ai],x[aj],dx);                   /*   3          */
    dr2        = iprod(dx,dx);                            /*   5          */
    
    if (dr2 == 0.0)
      continue;
      
    dr         = sqrt(dr2);                               /*  10          */
    dist       = dr-b0;
    kdist      = kb*dist;
    kdist2     = kdist*dist;
    
    vbond      = kdist2 + kcub*kdist2*dist;
    fbond      = -(two*kdist + three*kdist2*kcub)/dr;

    vtot      += vbond;       /* 21 */
    
    ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
    ki=IVEC2IS(dt);
    for (m=0; (m<DIM); m++) {                          /*  15          */
      fij=fbond*dx[m];
      f[ai][m]+=fij;
      f[aj][m]-=fij;
      fr->fshift[ki][m]+=fij;
      fr->fshift[CENTRAL][m]-=fij;
    }
  }                                           /*  54 TOTAL    */
  return vtot;
}

real harmonic(real kA,real kB,real xA,real xB,real x,real lambda,
	      real *V,real *F)
{
  const real half=0.5;
  real  L1,kk,x0,dx,dx2;
  real  v,f,dvdl;
  
  L1    = 1.0-lambda;
  kk    = L1*kA+lambda*kB;
  x0    = L1*xA+lambda*xB;
  
  dx    = x-x0;
  dx2   = dx*dx;
  
  f     = -kk*dx;
  v     = half*kk*dx2;
  dvdl  = half*(kB-kA)*dx2 + (xA-xB)*kk*dx;
  
  *F    = f;
  *V    = v;
  
  return dvdl;
  
  /* That was 19 flops */
}


real bonds(int nbonds,
	   t_iatom forceatoms[],t_iparams forceparams[],
	   rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	   matrix box,real lambda,real *dvdlambda,
	   t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  int  i,m,ki,ai,aj,type;
  real dr,dr2,fbond,vbond,fij,vtot;
  rvec dx;
  ivec dt;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
  
    pbc_rvec_sub(x[ai],x[aj],dx);			/*   3 		*/
    dr2=iprod(dx,dx);				/*   5		*/
    dr=sqrt(dr2);					/*  10		*/

    *dvdlambda += harmonic(forceparams[type].harmonic.krA,
			   forceparams[type].harmonic.krB,
			   forceparams[type].harmonic.rA,
			   forceparams[type].harmonic.rB,
			   dr,lambda,&vbond,&fbond);  /*  19  */

    if (dr2 == 0.0)
      continue;

    
    vtot  += vbond;/* 1*/
    fbond *= invsqrt(dr2);			/*   6		*/
#ifdef DEBUG
    if (debug)
      fprintf(debug,"BONDS: dr = %10g  vbond = %10g  fbond = %10g\n",
	      dr,vbond,fbond);
#endif
    ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
    ki=IVEC2IS(dt);

    for (m=0; (m<DIM); m++) {			/*  15		*/
      fij=fbond*dx[m];
      f[ai][m]+=fij;
      f[aj][m]-=fij;
      fr->fshift[ki][m]+=fij;
      fr->fshift[CENTRAL][m]-=fij;
    }
  }					/* 59 TOTAL	*/
  return vtot;
}

real water_pol(int nbonds,
	       t_iatom forceatoms[],t_iparams forceparams[],
	       rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	       matrix box,real lambda,real *dvdlambda,
	       t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  /* This routine implements anisotropic polarizibility for water, through
   * a shell connected to a dummy with spring constant that differ in the
   * three spatial dimensions in the molecular frame.
   */
  int  i,m,aO,aH1,aH2,aD,aS,type;
  rvec dOH1,dOH2,dHH,dOD,dDS,nW,kk,dx,kdx,proj;
#ifdef DEBUG
  rvec df;
#endif
  real vtot,fij,r_HH,r_OD,r_nW,tx,ty,tz;
  
  vtot = 0.0;
  if (nbonds > 0) {
    type   = forceatoms[0];
    kk[XX] = forceparams[type].wpol.kx;
    kk[YY] = forceparams[type].wpol.ky;
    kk[ZZ] = forceparams[type].wpol.kz;
    r_HH   = 1.0/forceparams[type].wpol.rHH;
    r_OD   = 1.0/forceparams[type].wpol.rOD;
    if (debug) {
      fprintf(debug,"WPOL: kk  = %10.3f        %10.3f        %10.3f\n",
	      kk[XX],kk[YY],kk[ZZ]);
      fprintf(debug,"WPOL: rOH = %10.3f  rHH = %10.3f  rOD = %10.3f\n",
	      forceparams[type].wpol.rOH,
	      forceparams[type].wpol.rHH,
	      forceparams[type].wpol.rOD);
    }
    for(i=0; (i<nbonds); ) {
      type = forceatoms[i++];
      aO   = forceatoms[i++];
      aH1  = aO+1;
      aH2  = aO+2;
      aD   = aO+3;
      aS   = aO+4;
      
      /* Compute vectors describing the water frame */
      rvec_sub(x[aH1],x[aO], dOH1);
      rvec_sub(x[aH2],x[aO], dOH2);
      rvec_sub(x[aH2],x[aH1],dHH);
      rvec_sub(x[aD], x[aO], dOD);
      rvec_sub(x[aS], x[aD], dDS);
      oprod(dOH1,dOH2,nW);
      
      /* Compute inverse length of normal vector 
       * (this one could be precomputed, but I'm too lazy now)
       */
      r_nW = invsqrt(iprod(nW,nW));
      /* This is for precision, but does not make a big difference,
       * it can go later.
       */
      r_OD = invsqrt(iprod(dOD,dOD)); 
      
      /* Normalize the vectors in the water frame */
      svmul(r_nW,nW,nW);
      svmul(r_HH,dHH,dHH);
      svmul(r_OD,dOD,dOD);
      
      /* Compute displacement of shell along components of the vector */
      dx[ZZ] = iprod(dDS,dOD);
      /* Compute projection on the XY plane: dDS - dx[ZZ]*dOD */
      for(m=0; (m<DIM); m++)
	proj[m] = dDS[m]-dx[ZZ]*dOD[m];
      
      /*dx[XX] = iprod(dDS,nW);
	dx[YY] = iprod(dDS,dHH);*/
      dx[XX] = iprod(proj,nW);
      for(m=0; (m<DIM); m++)
	proj[m] -= dx[XX]*nW[m];
      dx[YY] = iprod(proj,dHH);
      /*#define DEBUG*/
#ifdef DEBUG
      if (debug) {
	fprintf(debug,"WPOL: dx2=%10g  dy2=%10g  dz2=%10g  sum=%10g  dDS^2=%10g\n",
		sqr(dx[XX]),sqr(dx[YY]),sqr(dx[ZZ]),iprod(dx,dx),iprod(dDS,dDS));
	fprintf(debug,"WPOL: dHH=(%10g,%10g,%10g)\n",dHH[XX],dHH[YY],dHH[ZZ]);
	fprintf(debug,"WPOL: dOD=(%10g,%10g,%10g), 1/r_OD = %10g\n",
		dOD[XX],dOD[YY],dOD[ZZ],1/r_OD);
	fprintf(debug,"WPOL: nW =(%10g,%10g,%10g), 1/r_nW = %10g\n",
		nW[XX],nW[YY],nW[ZZ],1/r_nW);
	fprintf(debug,"WPOL: dx  =%10g, dy  =%10g, dz  =%10g\n",
		dx[XX],dx[YY],dx[ZZ]);
 	fprintf(debug,"WPOL: dDSx=%10g, dDSy=%10g, dDSz=%10g\n",
		dDS[XX],dDS[YY],dDS[ZZ]);
      }
#endif
      /* Now compute the forces and energy */
      kdx[XX] = kk[XX]*dx[XX];
      kdx[YY] = kk[YY]*dx[YY];
      kdx[ZZ] = kk[ZZ]*dx[ZZ];
      vtot   += iprod(dx,kdx);
      for(m=0; (m<DIM); m++) {
	/* This is a tensor operation but written out for speed */
	tx        =  nW[m]*kdx[XX];
	ty        = dHH[m]*kdx[YY];
	tz        = dOD[m]*kdx[ZZ];
	fij       = -tx-ty-tz;
#ifdef DEBUG
	df[m] = fij;
#endif
	f[aS][m] += fij;
	f[aD][m] -= fij;
      }
#ifdef DEBUG
      if (debug) {
	fprintf(debug,"WPOL: vwpol=%g\n",0.5*iprod(dx,kdx));
	fprintf(debug,"WPOL: df = (%10g, %10g, %10g)\n",df[XX],df[YY],df[ZZ]);
      }
#endif
    }	
  }
  return 0.5*vtot;
}

real bond_angle(matrix box,
		rvec xi,rvec xj,rvec xk,	/* in  */
		rvec r_ij,rvec r_kj,real *costh)		/* out */
/* Return value is the angle between the bonds i-j and j-k */
{
  /* 41 FLOPS */
  real th;
  
  pbc_rvec_sub(xi,xj,r_ij);			/*  3		*/
  pbc_rvec_sub(xk,xj,r_kj);			/*  3		*/

  *costh=cos_angle(r_ij,r_kj);		/* 25		*/
  th=acos(*costh);			/* 10		*/
					/* 41 TOTAL	*/
  return th;
}

real angles(int nbonds,
	    t_iatom forceatoms[],t_iparams forceparams[],
	    rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	    matrix box,real lambda,real *dvdlambda,
	    t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  int  i,ai,aj,ak,t1,t2,type;
  rvec r_ij,r_kj;
  real cos_theta,theta,dVdt,va,vtot;
  ivec jt,dt_ij,dt_kj;
  
  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    
    theta  = bond_angle(box,x[ai],x[aj],x[ak],
			r_ij,r_kj,&cos_theta);	/*  41		*/
  
    *dvdlambda += harmonic(forceparams[type].harmonic.krA,
			   forceparams[type].harmonic.krB,
			   forceparams[type].harmonic.rA*DEG2RAD,
			   forceparams[type].harmonic.rB*DEG2RAD,
			   theta,lambda,&va,&dVdt);  /*  21  */
    vtot += va;
    
    {
      int  m;
      real snt,st,sth;
      real cik,cii,ckk;
      real nrkj2,nrij2;
      rvec f_i,f_j,f_k;
      
      snt=sin(theta);				/*  10		*/
      if (fabs(snt) < 1e-12)
	snt=1e-12;
      st  = dVdt/snt;				/*  10		*/
      sth = st*cos_theta;			/*   1		*/
#ifdef DEBUG
      if (debug)
	fprintf(debug,"ANGLES: theta = %10g  vth = %10g  dV/dtheta = %10g\n",
		theta*RAD2DEG,va,dVdt);
#endif
      nrkj2=iprod(r_kj,r_kj);			/*   5		*/
      nrij2=iprod(r_ij,r_ij);
      
      cik=st*invsqrt(nrkj2*nrij2);		/*  12		*/ 
      cii=sth/nrij2;				/*  10		*/
      ckk=sth/nrkj2;				/*  10		*/
      
      for (m=0; (m<DIM); m++) {			/*  39		*/
	f_i[m]=-(cik*r_kj[m]-cii*r_ij[m]);
	f_k[m]=-(cik*r_ij[m]-ckk*r_kj[m]);
	f_j[m]=-f_i[m]-f_k[m];
	f[ai][m]+=f_i[m];
	f[aj][m]+=f_j[m];
	f[ak][m]+=f_k[m];
      }
      copy_ivec(SHIFT_IVEC(g,aj),jt);
      
    ivec_sub(SHIFT_IVEC(g,ai),jt,dt_ij);
    ivec_sub(SHIFT_IVEC(g,ak),jt,dt_kj);
    t1=IVEC2IS(dt_ij);
    t2=IVEC2IS(dt_kj);
      
      rvec_inc(fr->fshift[t1],f_i);
      rvec_inc(fr->fshift[CENTRAL],f_j);
      rvec_inc(fr->fshift[t2],f_k);
    }                                           /* 168 TOTAL	*/
  }
  return vtot;
}

real dih_angle(matrix box,
	       rvec xi,rvec xj,rvec xk,rvec xl,
	       rvec r_ij,rvec r_kj,rvec r_kl,rvec m,rvec n,
	       real *cos_phi,real *sign)
{
  real ipr,phi;

  pbc_rvec_sub(xi,xj,r_ij);       		/*  3 		*/
  pbc_rvec_sub(xk,xj,r_kj);			/*  3		*/
  pbc_rvec_sub(xk,xl,r_kl);			/*  3		*/

  oprod(r_ij,r_kj,m); 			/*  9 		*/
  oprod(r_kj,r_kl,n);			/*  9		*/
  *cos_phi=cos_angle(m,n); 		/* 41 		*/
  phi=acos(*cos_phi); 			/* 10 		*/
  ipr=iprod(r_ij,n); 			/*  5 		*/
  (*sign)=(ipr<0.0)?-1.0:1.0;
  phi=(*sign)*phi; 			/*  1		*/
					/* 84 TOTAL	*/
  return phi;
}

void do_dih_fup(int i,int j,int k,int l,real ddphi,
		rvec r_ij,rvec r_kj,rvec r_kl,
		rvec m,rvec n,rvec f[],t_forcerec *fr,t_graph *g,
		rvec x[])
{
  /* 143 FLOPS */
  rvec f_i,f_j,f_k,f_l;
  rvec uvec,vvec,svec;
  real ipr,nrkj,nrkj2;
  real a,p,q;
  int  t1,t2,t3;
  ivec jt,dt_ij,dt_kj,dt_lj;  
  
  ipr   = iprod(m,m);		/*  5 	*/
  nrkj2 = iprod(r_kj,r_kj);	/*  5	*/
  nrkj  = sqrt(nrkj2);		/* 10	*/
  a     = -ddphi*nrkj/ipr;	/* 11	*/
  svmul(a,m,f_i);		/*  3	*/
  ipr   = iprod(n,n);		/*  5	*/
  a     = ddphi*nrkj/ipr;	/* 11	*/
  svmul(a,n,f_l);		/*  3 	*/
  p     = iprod(r_ij,r_kj);	/*  5	*/
  p    /= nrkj2;		/* 10	*/
  q     = iprod(r_kl,r_kj);	/*  5	*/
  q    /= nrkj2;		/* 10	*/
  svmul(p,f_i,uvec);		/*  3	*/
  svmul(q,f_l,vvec);		/*  3	*/
  rvec_sub(uvec,vvec,svec);	/*  3	*/
  rvec_sub(f_i,svec,f_j);	/*  3	*/
  rvec_add(f_l,svec,f_k);	/*  3	*/
  rvec_inc(f[i],f_i);   	/*  3	*/
  rvec_dec(f[j],f_j);   	/*  3	*/
  rvec_dec(f[k],f_k);   	/*  3	*/
  rvec_inc(f[l],f_l);   	/*  3	*/

  
  copy_ivec(SHIFT_IVEC(g,j),jt);
  ivec_sub(SHIFT_IVEC(g,i),jt,dt_ij);
  ivec_sub(SHIFT_IVEC(g,k),jt,dt_kj);
  ivec_sub(SHIFT_IVEC(g,l),jt,dt_lj);
  t1=IVEC2IS(dt_ij);
  t2=IVEC2IS(dt_kj);
  t3=IVEC2IS(dt_lj);
  
  rvec_inc(fr->fshift[t1],f_i);
  rvec_dec(fr->fshift[CENTRAL],f_j);
  rvec_dec(fr->fshift[t2],f_k);
  rvec_inc(fr->fshift[t3],f_l);
  /* 112 TOTAL 	*/
}


real dopdihs(real cpA,real cpB,real phiA,real phiB,int mult,
	     real phi,real lambda,real *V,real *F)
{
  real v,dvdl,mdphi,v1,sdphi,ddphi;
  real L1   = 1.0-lambda;
  real ph0  = DEG2RAD*(L1*phiA+lambda*phiB);
  real cp   = L1*cpA + lambda*cpB;
  
  mdphi =  mult*phi-ph0;
  sdphi = sin(mdphi);
  ddphi = -cp*mult*sdphi;
  v1    = 1.0+cos(mdphi);
  v     = cp*v1;
  
  dvdl  = (cpB-cpA)*v1 - cp*(phiA-phiB)*sdphi;
  
  *V = v;
  *F = ddphi;
  
  return dvdl;
  
  /* That was 40 flops */
}

static real dopdihs_min(real cpA,real cpB,real phiA,real phiB,int mult,
		       real phi,real lambda,real *V,real *F)
     /* similar to dopdihs, except for a minus sign  *
      * and a different treatment of mult/phi0       */
{
  real v,dvdl,mdphi,v1,sdphi,ddphi;
  real L1   = 1.0-lambda;
  real ph0  = DEG2RAD*(L1*phiA+lambda*phiB);
  real cp   = L1*cpA + lambda*cpB;
  
  mdphi = mult*(phi-ph0);
  sdphi = sin(mdphi);
  ddphi = cp*mult*sdphi;
  v1    = 1.0-cos(mdphi);
  v     = cp*v1;
  
  dvdl  = (cpB-cpA)*v1 - cp*(phiA-phiB)*sdphi;
  
  *V = v;
  *F = ddphi;
  
  return dvdl;
  
  /* That was 40 flops */
}

real pdihs(int nbonds,
	   t_iatom forceatoms[],t_iparams forceparams[],
	   rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	   matrix box,real lambda,real *dvdlambda,
	   t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  int  i,type,ai,aj,ak,al;
  rvec r_ij,r_kj,r_kl,m,n;
  real phi,cos_phi,sign,ddphi,vpd,vtot;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    al   = forceatoms[i++];
    
    phi=dih_angle(box,x[ai],x[aj],x[ak],x[al],r_ij,r_kj,r_kl,m,n,
		  &cos_phi,&sign);			/*  84 		*/
		
    *dvdlambda += dopdihs(forceparams[type].pdihs.cpA,
			  forceparams[type].pdihs.cpB,
			  forceparams[type].pdihs.phiA,
			  forceparams[type].pdihs.phiB,
			  forceparams[type].pdihs.mult,
			  phi,lambda,&vpd,&ddphi);
		       
    vtot += vpd;
    do_dih_fup(ai,aj,ak,al,ddphi,r_ij,r_kj,r_kl,m,n,
	       f,fr,g,x); 				/* 112		*/

#ifdef DEBUG
    fprintf(debug,"pdih: (%d,%d,%d,%d) cp=%g, phi=%g\n",
	    ai,aj,ak,al,cos_phi,phi);
#endif
  } /* 223 TOTAL 	*/

  return vtot;
}

real idihs(int nbonds,
	   t_iatom forceatoms[],t_iparams forceparams[],
	   rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	   matrix box,real lambda,real *dvdlambda,
	   t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  int  i,type,ai,aj,ak,al;
  real phi,cos_phi,ddphi,sign,vid,vtot;
  rvec r_ij,r_kj,r_kl,m,n;
  
  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    al   = forceatoms[i++];
    
    phi=dih_angle(box,x[ai],x[aj],x[ak],x[al],r_ij,r_kj,r_kl,m,n,
		  &cos_phi,&sign);			/*  84		*/
    
    *dvdlambda += harmonic(forceparams[type].harmonic.krA,
			   forceparams[type].harmonic.krB,
			   forceparams[type].harmonic.rA*DEG2RAD,
			   forceparams[type].harmonic.rB*DEG2RAD,
			   phi,lambda,&vid,&ddphi);    /*   21          */   

    vtot += vid;
    do_dih_fup(ai,aj,ak,al,(real)(-ddphi),r_ij,r_kj,r_kl,m,n,
	       f,fr,g,x);				/* 112		*/
    /* 217 TOTAL	*/
#ifdef DEBUG
    fprintf("idih: (%d,%d,%d,%d) cp=%g, phi=%g\n",
	    ai,aj,ak,al,cos_phi,phi);
#endif
  }
  return vtot;
}

real posres(int nbonds,
	    t_iatom forceatoms[],t_iparams forceparams[],
	    rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	    matrix box,real lambda,real *dvdlambda,
	    t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  int  i,ai,m,type;
  real v,vtot,fi,*fc;
  rvec dx;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    fc   = forceparams[type].posres.fc;

    if (fr->ePBC == epbcNONE)
      rvec_sub(x[ai],forceparams[type].posres.pos0,dx);
    else
      pbc_dx(x[ai],forceparams[type].posres.pos0,dx);
    v=0;
    for (m=0; (m<DIM); m++) {
      fi        = f[ai][m] - fc[m]*dx[m];
      v        += 0.5*fc[m]*dx[m]*dx[m];
      f[ai][m]  = fi;
    }
    vtot += v;
  }
  return vtot;
}

static real low_angres(int nbonds,
		       t_iatom forceatoms[],t_iparams forceparams[],
		       rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
		       matrix box,real lambda,real *dvdlambda,
		       bool bZAxis)
{
  int  i,m,type,ai,aj,ak,al,t;
  real phi,cos_phi,vid,vtot,dVdphi;
  rvec r_ij,r_kl,f_i,f_k;
  real st,sth,sin_phi,nrij2,nrkl2,c,cij,ckl;

  ivec it,jt,dt;  
  
  vtot = 0.0;
  ak=al=0; /* to avoid warnings */
  for(i=0; i<nbonds; ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    pbc_rvec_sub(x[aj],x[ai],r_ij);            	/*  3		*/
    if (!bZAxis) {      
      ak   = forceatoms[i++];
      al   = forceatoms[i++];
      pbc_rvec_sub(x[al],x[ak],r_kl);           /*  3		*/
    } else {
      r_kl[XX] = 0;
      r_kl[YY] = 0;
      r_kl[ZZ] = 1;
    }

    cos_phi = cos_angle(r_ij,r_kl);		/* 25		*/
    phi     = acos(cos_phi);                    /* 10           */

    *dvdlambda += dopdihs_min(forceparams[type].pdihs.cpA,
			      forceparams[type].pdihs.cpB,
			      forceparams[type].pdihs.phiA,
			      forceparams[type].pdihs.phiB,
			      forceparams[type].pdihs.mult,
			      phi,lambda,&vid,&dVdphi); /*  40  */
    
    vtot += vid;

    sin_phi = sin(phi);			        /*  10		*/
    if (fabs(sin_phi) < 1e-12)
      sin_phi=1e-12;
    st  = -dVdphi/sin_phi;		        /*  10		*/
    sth = st*cos_phi;				/*   1		*/
    nrij2 = iprod(r_ij,r_ij);			/*   5		*/
    nrkl2 = iprod(r_kl,r_kl);                   /*   5          */
    
    c   = st*invsqrt(nrij2*nrkl2);		/*  11		*/ 
    cij = sth/nrij2;				/*  10		*/
    ckl = sth/nrkl2;				/*  10		*/
    
    for (m=0; m<DIM; m++) {			/*  18+18       */
      f_i[m] = (c*r_kl[m]-cij*r_ij[m]);
      f[ai][m] += f_i[m];
      f[aj][m] -= f_i[m];
      if (!bZAxis) {
	f_k[m] = (c*r_ij[m]-ckl*r_kl[m]);
	f[ak][m] += f_k[m];
	f[al][m] -= f_k[m];
      }
    }
    
    ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
    t=IVEC2IS(dt);
    
    rvec_inc(fr->fshift[t],f_i);
    rvec_dec(fr->fshift[CENTRAL],f_i);
    if (!bZAxis) {
    ivec_sub(SHIFT_IVEC(g,ak),SHIFT_IVEC(g,al),dt);
    ivec_sub(it,jt,dt);
    t=IVEC2IS(dt);
      rvec_inc(fr->fshift[t],f_k);
      rvec_dec(fr->fshift[CENTRAL],f_k);
    }
  }

  return vtot;  /*  191 / 164 (bZAxis)  total  */
}

real angres(int nbonds,
	    t_iatom forceatoms[],t_iparams forceparams[],
	    rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	    matrix box,real lambda,real *dvdlambda,
	    t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  return low_angres(nbonds,forceatoms,forceparams,x,f,fr,g,box,
		    lambda,dvdlambda,FALSE);
}

real angresz(int nbonds,
	     t_iatom forceatoms[],t_iparams forceparams[],
	     rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	     matrix box,real lambda,real *dvdlambda,
	     t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  return low_angres(nbonds,forceatoms,forceparams,x,f,fr,g,box,
		    lambda,dvdlambda,TRUE);
}

real unimplemented(int nbonds,
		   t_iatom forceatoms[],t_iparams forceparams[],
		   rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
		   matrix box,real lambda,real *dvdlambda,
		   t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  fatal_error(0,"*** you are using a not implemented function");

  return 0.0; /* To make the compiler happy */
}

static void my_fatal(char *s,int ai,int aj,int ak,int al)
{
  fatal_error(0,"%s is NaN in rbdih, ai-ak=%d,%d,%d,%d",s,ai,aj,ak,al);
}

#define CHECK(s) if ((s != s) || (s < -1e10) || (s > 1e10)) my_fatal(#s,ai,aj,ak,al)

real rbdihs(int nbonds,
	    t_iatom forceatoms[],t_iparams forceparams[],
	    rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	    matrix box,real lambda,real *dvdlambda,
	    t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  static const real c0=0.0,c1=1.0,c2=2.0,c3=3.0,c4=4.0,c5=5.0;
  int  type,ai,aj,ak,al,i,j;
  rvec r_ij,r_kj,r_kl,m,n;
  real parm[NR_RBDIHS];
  real phi,cos_phi,rbp;
  real v,sign,ddphi,sin_phi;
  real cosfac,vtot;

  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    al   = forceatoms[i++];

    phi=dih_angle(box,x[ai],x[aj],x[ak],x[al],r_ij,r_kj,r_kl,m,n,
		  &cos_phi,&sign);			/*  84		*/

    /* Change to polymer convention */
    if (phi < c0)
      phi += M_PI;
    else
      phi -= M_PI;			/*   1		*/
    cos_phi = -cos_phi;					/*   1		*/
    
    sin_phi=sin(phi);
    
    for(j=0; (j<NR_RBDIHS); j++)
      parm[j] = forceparams[type].rbdihs.rbc[j];
    
    /* Calculate cosine powers */
    /* Calculate the energy */
    /* Calculate the derivative */
    v       = parm[0];
    ddphi   = c0;
    cosfac  = c1;
    
    rbp     = parm[1];
    ddphi  += rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    rbp     = parm[2];
    ddphi  += c2*rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    rbp     = parm[3];
    ddphi  += c3*rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    rbp     = parm[4];
    ddphi  += c4*rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    rbp     = parm[5];
    ddphi  += c5*rbp*cosfac;
    cosfac *= cos_phi;
    v      += cosfac*rbp;
    
    ddphi = -ddphi*sin_phi;				/*  11		*/
    
    do_dih_fup(ai,aj,ak,al,ddphi,r_ij,r_kj,r_kl,m,n,
	       f,fr,g,x);				/* 112		*/
    vtot += v;
  }  
  return vtot;
}

/***********************************************************
 *
 *   G R O M O S  9 6   F U N C T I O N S
 *
 ***********************************************************/
real g96harmonic(real kA,real kB,real xA,real xB,real x,real lambda,
		 real *V,real *F)
{
  const real half=0.5;
  real  L1,kk,x0,dx,dx2;
  real  v,f,dvdl;
  
  L1    = 1.0-lambda;
  kk    = L1*kA+lambda*kB;
  x0    = L1*xA+lambda*xB;
  
  dx    = x-x0;
  dx2   = dx*dx;
  
  f     = -kk*dx;
  v     = half*kk*dx2;
  dvdl  = half*(kB-kA)*dx2 + (xA-xB)*kk*dx;
  
  *F    = f;
  *V    = v;
  
  return dvdl;
  
  /* That was 21 flops */
}

real g96bonds(int nbonds,
	      t_iatom forceatoms[],t_iparams forceparams[],
	      rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	      matrix box,real lambda,real *dvdlambda,
	      t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  int  i,m,ki,ai,aj,type;
  real dr2,fbond,vbond,fij,vtot;
  rvec dx;
  ivec dt;
  
  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
  
    pbc_rvec_sub(x[ai],x[aj],dx);		/*   3 		*/
    dr2=iprod(dx,dx);				/*   5		*/
      
    *dvdlambda += g96harmonic(forceparams[type].harmonic.krA,
			      forceparams[type].harmonic.krB,
			      forceparams[type].harmonic.rA,
			      forceparams[type].harmonic.rB,
			      dr2,lambda,&vbond,&fbond);

    vtot  += 0.5*vbond;                             /* 1*/
#ifdef DEBUG
    if (debug)
      fprintf(debug,"G96-BONDS: dr = %10g  vbond = %10g  fbond = %10g\n",
	      sqrt(dr2),vbond,fbond);
#endif
   
    ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);
    ki=IVEC2IS(dt);
    for (m=0; (m<DIM); m++) {			/*  15		*/
      fij=fbond*dx[m];
      f[ai][m]+=fij;
      f[aj][m]-=fij;
      fr->fshift[ki][m]+=fij;
      fr->fshift[CENTRAL][m]-=fij;
    }
  }					/* 44 TOTAL	*/
  return vtot;
}

real g96bond_angle(matrix box,
		   rvec xi,rvec xj,rvec xk,	/* in  */
		   rvec r_ij,rvec r_kj)		/* out */
/* Return value is the angle between the bonds i-j and j-k */
{
  real costh;
  
  pbc_rvec_sub(xi,xj,r_ij);			/*  3		*/
  pbc_rvec_sub(xk,xj,r_kj);			/*  3		*/

  costh=cos_angle(r_ij,r_kj);		/* 25		*/
					/* 41 TOTAL	*/
  return costh;
}

real g96angles(int nbonds,
	       t_iatom forceatoms[],t_iparams forceparams[],
	       rvec x[],rvec f[],t_forcerec *fr,t_graph *g,
	       matrix box,real lambda,real *dvdlambda,
	       t_mdatoms *md,int ngrp,real egnb[],real egcoul[])
{
  int  i,ai,aj,ak,type,m,t1,t2;
  rvec r_ij,r_kj;
  real cos_theta,dVdt,va,vtot;
  real rij_1,rij_2,rkj_1,rkj_2,rijrkj_1;
  rvec f_i,f_j,f_k;
  ivec jt,dt_ij,dt_kj;
  
  vtot = 0.0;
  for(i=0; (i<nbonds); ) {
    type = forceatoms[i++];
    ai   = forceatoms[i++];
    aj   = forceatoms[i++];
    ak   = forceatoms[i++];
    
    cos_theta  = g96bond_angle(box,x[ai],x[aj],x[ak],r_ij,r_kj);	
    
    *dvdlambda += g96harmonic(forceparams[type].harmonic.krA,
			      forceparams[type].harmonic.krB,
			      forceparams[type].harmonic.rA,
			      forceparams[type].harmonic.rB,
			      cos_theta,lambda,&va,&dVdt);
    vtot    += va;
    
    rij_1    = invsqrt(iprod(r_ij,r_ij));
    rkj_1    = invsqrt(iprod(r_kj,r_kj));
    rij_2    = rij_1*rij_1;
    rkj_2    = rkj_1*rkj_1;
    rijrkj_1 = rij_1*rkj_1;                     /* 23 */
    
#ifdef DEBUG
    if (debug)
      fprintf(debug,"G96ANGLES: costheta = %10g  vth = %10g  dV/dct = %10g\n",
	      cos_theta,va,dVdt);
#endif
    for (m=0; (m<DIM); m++) {			/*  42	*/
      f_i[m]=dVdt*(r_kj[m]*rijrkj_1 - r_ij[m]*rij_2*cos_theta);
      f_k[m]=dVdt*(r_ij[m]*rijrkj_1 - r_kj[m]*rkj_2*cos_theta);
      f_j[m]=-f_i[m]-f_k[m];
      f[ai][m]+=f_i[m];
      f[aj][m]+=f_j[m];
      f[ak][m]+=f_k[m];
    }
    copy_ivec(SHIFT_IVEC(g,aj),jt);
    
    ivec_sub(SHIFT_IVEC(g,ai),jt,dt_ij);
    ivec_sub(SHIFT_IVEC(g,ak),jt,dt_kj);
    t1=IVEC2IS(dt_ij);
    t2=IVEC2IS(dt_kj);

    rvec_inc(fr->fshift[t1],f_i);
    rvec_inc(fr->fshift[CENTRAL],f_j);
    rvec_inc(fr->fshift[t2],f_k);               /* 9 */
    /* 163 TOTAL	*/
  }
  return vtot;
}

