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
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_ewald_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "vec.h"
#include "complex.h"
#include "smalloc.h"
#include "futil.h"
#include "lrutil.h"
#include "fftgrid.h"
#include "pppm.h"
#include "fatal.h"
/*#include "ghat.h"*/

#define TOL 2e-5

#define NCOS 360
static real ncos[NCOS+1];
static real nsin[NCOS+1];
static real twopi  = 2.0*M_PI;
static real two_pi = 0.5/M_PI;

real mysincos(real x,real sincos[])
{
  real ix,dx,xtp;
  int  i;
  
  /*while(x < 0)
    x+=twopi;
  while(x>twopi)
  x-=twopi;*/
  xtp = x*two_pi;
  x  = (xtp - (int) xtp);
  ix = (x*NCOS);
  i  = ix;
  dx = ix-i;
  
  return (sincos[i]*(1.0-dx) + dx*(sincos[i+1]));
}

#define mycos(x) mysincos(x,ncos)
#define mysin(x) mysincos(x,nsin)

static void init_sincos(void)
{
  int  i;
  real x;
  
  for(i=0; (i<=NCOS); i++) {
    x = (i*twopi)/NCOS; 
    ncos[i] = cos(x);
    nsin[i] = sin(x);
  }
}

static t_complex Qk(rvec k,int nj,rvec rj[],real qj[])
/* Return fourier component of discrete charge distribution for all particles
 * (Eq. 33)
 */
{
  int       i;
  real      ipr;
  t_complex c,cqi,cpk;
  
  c = cnul;
  for(i=0; (i<nj); i++) {
    if (qj[i] != 0.0) {
      ipr = iprod(k,rj[i]);
      cpk = rcexp(-ipr);
      cqi = rcmul(qj[i],cpk);
      c   = cadd(c,cqi);
    }
  }
  /*if (fabs(c.re) < TOL) 
    c.re = 0;
    if (fabs(c.im) < TOL)
    c.im = 0;*/
  
  return c;
}

static t_complex exp_ikr(t_complex eix,t_complex eiy,t_complex eiz)
{
  /* Function derived using maple, see complex.map */
  t_complex c;
  
  c.re = ((eix.re*(eiy.re*eiz.re-eiy.im*eiz.im)) -
	  (eix.im*(eiy.re*eiz.im+eiy.im*eiz.re)));
  c.im = ((eix.re*(eiy.re*eiz.im+eiy.im*eiz.re)) +
	  (eix.im*(eiy.re*eiz.re-eiy.im*eiz.im)));
  if (debug) 
    fprintf(debug,"c.re: %g, c.im: %g\n",c.re,c.im);
  
  return c;
}

static t_complex exp_min_ikr(t_complex eix,t_complex eiy,t_complex eiz)
{
  /* Function derived using maple, see complex.map */
  t_complex c;
  
  c.re = ((eix.re*eiy.re-eix.im*eiy.im)*eiz.re -
	  (eix.im*eiy.re+eix.re*eiy.im)*eiz.im);
  c.im = ((eix.im*eiy.im-eix.re*eiy.re)*eiz.im -
	  (eix.im*eiy.re+eix.re*eiy.im)*eiz.re);
  if (debug) 
    fprintf(debug,"c.re: %g, c.im: %g\n",c.re,c.im);
  
  return c;
}

static t_complex Qk2(int nj,rvec rj[],real qj[],cvec *eix,cvec *eiy,cvec *eiz)
/* Return fourier component of discrete charge distribution for all particles
 * (Eq. 33)
 */
{
  int       i;
  t_complex c;
  
  c = cnul;
  for(i=0; (i<nj); i++) 
    if (qj[i] != 0.0) 
      c = cadd(c,rcmul(qj[i],exp_min_ikr(eix[i][XX],eiy[i][YY],eiz[i][ZZ])));
  
  return c;
}

static void mk_qtab(int nx,int ny,int nz,t_complex ***qtab,
		    rvec box,int nj,rvec rj[],real charge[])
{
  int       ix,iy,iz;
  int       m;
  rvec      k,lll;
  t_complex c;
  
  calc_lll(box,lll);
  
  /* Loop over lattice vectors in fourier space */    
  for(ix=0; (ix < nx); ix++) {
    for(iy=0; (iy < ny); iy++) {
      for(iz=0; (iz < nz); iz++) {
	if ((ix == 0) && (iy == 0) && (iz == 0)) 
	  c = cnul;
	else {
	  calc_k(lll,ix,iy,iz,nx,ny,nz,k);
	  
	  c = Qk(k,nj,rj,charge);
	}
	qtab[ix][iy][iz] = c;
      }
    }
  }
}

static t_complex phi_k(rvec k,rvec rj,real gk,t_complex qk)
{
  real      phikr;
  t_complex cp,cpkr;
  
  phikr  = iprod(k,rj);
  cpkr   = rcexp(phikr);
  cp     = rcmul(gk,cmul(qk,cpkr));
  
  return cp;
}

static t_complex phi_k2(real gk,t_complex qk,cvec eix,cvec eiy,cvec eiz)
{
  return rcmul(gk,cmul(qk,exp_ikr(eix[XX],eiy[YY],eiz[ZZ])));
}

void print_qkgrid(FILE *out,char *title,t_complex ***grid,real factor,
		  char *pdb,rvec box,int nx,int ny,int nz)
{
  static char *pdbformat="%-6s%5d  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
  FILE     *fp;
  int      i,ix,iy,iz;
  real     fac=50.0;
  rvec     boxfac;
  t_complex g;
  
  if (pdb)
    fp = ffopen(pdb,"w");
  else
    fp = out;
  if (!fp)
    return;

  boxfac[XX] = fac*box[XX]/nx;
  boxfac[YY] = fac*box[YY]/ny;
  boxfac[ZZ] = fac*box[ZZ]/nz;
  
  fprintf(fp,"Printing all non-zero complex elements of %s\n",title);
  for(i=ix=0; (ix<nx); ix++)
    for(iy=0; (iy<ny); iy++)
      for(iz=0; (iz<nz); iz++,i++) {
	g = grid[ix][iy][iz];
	if (fabs(g.im) > TOL) {
	  if (pdb) {
	    fprintf(fp,pdbformat,"ATOM",i,"H","H",' ',
		    i,ix*boxfac[XX],iy*boxfac[YY],iz*boxfac[ZZ],
		    1.0,g.im);
	  } 
	  else {
	    fprintf(fp,"%s[%2d][%2d][%2d] = %12.5e + i %12.5e%s\n",
		    title,ix,iy,iz,g.re*factor,g.im*factor,
		    (g.im != 0) ? " XXX" : "");
	  }
	}
      }
  fflush(fp);
}

static int mk_ivecs(int nx,int ny,int nz,ivec **ivptr)
{
  int  i,j,k,niv;
  ivec *iv;
  
  niv = 0;
  snew(iv,(2*nx+1)*(2*ny+1)*(2*nz+1));
  for(i=-nx; (i<=nx); i++)
    for(j=-ny; (j<=ny); j++)
      for(k=-nz; (k<=nz); k++) 
	if ((i!=0) || (j!=0) || (k!=0)) { 
	  iv[niv][XX] = i;
	  iv[niv][YY] = j;
	  iv[niv][ZZ] = k;
	  niv++;
	}
  *ivptr = iv;
  
  return niv;
}

static rvec box_save;

real box_length(ivec iv,rvec box)
{
  int m;
  rvec x;
  
  for(m=0; (m<DIM); m++) 
    x[m] = box[m]*(abs(iv[m]));
  return iprod(x,x);
}

static int box_cmp(const void *a,const void *b)
{
  int  *ia,*ib;
  real xa,xb;
  int  m;
  
  ia = (int *)a;
  ib = (int *)b;
  xa = box_length(ia,box_save);
  xb = box_length(ib,box_save);
  
  if (xa < xb)
    return -1;
  else if (xa > xb)
    return 1;
  else
    return 0;
}

static int sort_ivecs(int niv,ivec iv[],rvec box,int index[])
{
  const real tolerance = 1e-3;
  int   i,nind;
  real  length[2];
  int   cur=0;
#define new (1-cur)

  /* First step, sort all the boxes on distance from the origin */
  copy_rvec(box,box_save);
  qsort(iv,niv,sizeof(iv[0]),box_cmp);
  
  /* Next step sort all the indices so that we add symmetrically all boxes
   * at the same distance from the origin
   */
  nind        = 0;
  index[0]    = 0;
  length[cur] = box_length(iv[0],box);
  for(i=1; (i<niv); i++) {
    length[new] = box_length(iv[i],box);
    if ((length[new]/length[cur]) - 1 > tolerance) {
      index[++nind] = i;
      cur = new;
    }
  }
  index[++nind] = i;
  
  return nind;
}

static void print_ivecs(FILE *fp,int niv,ivec iv[],rvec box,
			int nind,int index[])
{
  int i,j;
 
  fprintf(fp,"EWALD BOX ORDER\n"); 
  fprintf(fp,"There are %d ivecs and %d indices\n",niv,nind);
  fprintf(fp,"Index   Ivec      X      Y      Z      Dist\n");
  for(i=0; (i<nind); i++) 
    for(j=index[i]; (j<index[i+1]); j++) {
      fprintf(fp,"%5d  %5d  %5d  %5d  %5d  %8.3f\n",
	      i,j,iv[j][XX],iv[j][YY],iv[j][ZZ],box_length(iv[j],box));
    }
}

static void mk_kprops(int niv,ivec iv[],int nx,int ny,int nz,rvec lll,
		      real rc,rvec all_k[],real gk[])
{
  int  i;
  
  for(i=0; (i<niv); i++) {
    calc_k(lll,iv[i][XX],iv[i][YY],iv[i][ZZ],nx,ny,nz,all_k[i]);
    gk[i] = gknew(norm(all_k[i]),rc,0.0);
  }
}

static void tabulate_eir(int natom,rvec x[],int kmax,cvec **eir,rvec lll)
{
  int  i,j,m;
  
  if (kmax < 1)
    fatal_error(0,"Go away! kmax = %d\n",kmax);
  for(i=0; (i<natom); i++) {
    for(m=0; (m<DIM); m++) {
      eir[0][i][m].re = 1;
      eir[0][i][m].im = 0;
    }
    for(m=0; (m<DIM); m++) {
      eir[1][i][m].re = cos(x[i][m]*lll[m]);
      eir[1][i][m].im = sin(x[i][m]*lll[m]); 
    }
    for(j=2; (j<kmax); j++) 
      for(m=0; (m<DIM); m++)
	eir[j][i][m] = cmul(eir[j-1][i][m],eir[1][i][m]);
  }
}

real do_ewald_new(FILE *log,       t_inputrec *ir,
		  int natoms,      rvec x[],rvec f[],
		  real charge[],   rvec box,
		  real phi[],      t_commrec *cr,
		  bool bOld)
{
  static    bool bFirst = TRUE;
  static    int       nx,ny,nz,kmax;
  static    ivec      *iv;
  static    int       *iv_index,nind,niv;
  static    rvec      *all_k;
  static    real      *gk;
  static    t_complex *qk;
  static    cvec      **eir;
  
  real      Gk,energy;
  int       m,ix,iy,iz,i,j,n;
  rvec      k,lll;
  t_complex phik,phitot,Ei[3],QQk;

  /* Now start the actual solution! */
  calc_lll(box,lll);
  
  if (bFirst) {  
    if (cr != NULL) {
      if (cr->nprocs > 1)
	fatal_error(0,"No parallel Ewald (yet)");
    }
    
    nx = ir->nkx;
    ny = ir->nky;
    nz = ir->nkz;
    
    fprintf(log,"Will do Ewald summation using %d x %d x %d k-vectors\n",
	    nx,ny,nz);
    fprintf(log,"This code is not optimized, so be prepared to wait!\n"); 
       
    niv   = mk_ivecs(nx/2-1,ny/2-1,nz/2-1,&iv);
    snew(iv_index,niv);
    nind  = sort_ivecs(niv,iv,box,iv_index);
    snew(all_k,niv);
    snew(gk,niv);
    mk_kprops(niv,iv,nx,ny,nz,lll,ir->rcoulomb,all_k,gk);
    snew(qk,niv);
    
    if (debug) 
      print_ivecs(debug,niv,iv,box,nind,iv_index);

    kmax = 1+max(nx/2,max(ny/2,nz/2));
    snew(eir,kmax);
    for(j=0; (j<kmax); j++)
      snew(eir[j],natoms);
            
    bFirst = FALSE;
  }
  
  /* Compute the exp(-ir) vectors */
  tabulate_eir(natoms,x,kmax,eir,lll);
  
  /* Compute the charge spread function */
  for(j=0; (j<nind); j++) {
    for(n=iv_index[j]; (n<iv_index[j+1]); n++) {
      qk[n] = Qk2(natoms,x,charge,
		  eir[iv[n][XX]],eir[iv[n][YY]],eir[iv[n][ZZ]]);
    }
  }
  
  energy = 0;
  for(i=0; (i<natoms); i++) {
    /* Loop over atoms */
    phitot = cnul;
    
    /* Electric field vector */
    for(m=0; (m<DIM); m++)
      Ei[m] = cnul;
      
    for(j=0; (j<nind); j++) {
      for(n=iv_index[j]; (n<iv_index[j+1]); n++) {
	QQk = qk[n];
	Gk  = gk[n];
	
	/* Calculate potential */
	phik   = phi_k2(Gk,QQk,eir[iv[n][XX]][i],eir[iv[n][YY]][i],
			eir[iv[n][ZZ]][i]);
	phitot = cadd(phik,phitot);
	
	for(m=0; (m<DIM); m++)
	  Ei[m] = cadd(Ei[m],rcmul(all_k[n][m],phik));
      }
    }
    /* Potential at atom i, and energy contribution */
#ifdef DEBUG
    fprintf(log,"phi[%3d] = %10.3e + i %10.3e\n",i,phitot.re,phitot.im);
#endif
    phi[i]  = phitot.re;
    energy += phi[i]*charge[i];
    
    /* Force on the atom */
    for(m=0; (m<DIM); m++) {
      f[i][m] = charge[i]*Ei[m].im;
#ifdef DEBUG
      if (fabs(Ei[m].re) > 1e-6)
	fprintf(log,"Ei[%1d] atom %3d = %10.3e + i%10.3e\n",
		m,i,Ei[m].re,Ei[m].im);
#endif
    }
  }
  return energy;
}


real do_ewald(FILE *log,       t_inputrec *ir,
	      int natoms,      rvec x[],rvec f[],
	      real charge[],   rvec box,
	      real phi[],      t_commrec *cr,
	      bool bOld)
{
  static    bool bFirst = TRUE;
  static    t_complex ***qtab;
  static    real      ***gtab;
  static    int       nx,ny,nz,kmax;
  static    ivec      *iv;
  static    int       *iv_index,nind,niv;
  static    rvec      *all_k;
  static    real      *gk;
  static    t_complex *qk;
  
  real      Gk,energy;
  int       m,ix,iy,iz,i,j,n;
  rvec      k,lll;
  t_complex phik,phitot,Ei[3],QQk;

  /* Now start the actual solution! */
  calc_lll(box,lll);
  
  if (bFirst) {  
    if (cr != NULL) {
      if (cr->nprocs > 1)
	fatal_error(0,"No parallel Ewald (yet)");
    }
    
    nx = ir->nkx;
    ny = ir->nky;
    nz = ir->nkz;
    
    fprintf(log,"Will do Ewald summation using %d x %d x %d k-vectors\n",
	    nx,ny,nz);
    fprintf(log,"This code is not optimized, so be prepared to wait!\n"); 
       
    gtab = mk_rgrid(nx,ny,nz);
    qtab = mk_cgrid(nx,ny,nz);

    mk_ghat(NULL,nx,ny,nz,gtab,box,ir->rshort,ir->rcoulomb,TRUE,bOld);
            
    bFirst = FALSE;
  }
  /* Tabulate charge spread & charge distribution in fourier space */
  mk_qtab(nx,ny,nz,qtab,box,natoms,x,charge);
  
  print_qkgrid(log,"qk",qtab,1.0,"qkfour.pdb",box,nx,ny,nz);
  
  energy = 0;
  for(i=0; (i<natoms); i++) {
    /* Loop over atoms */
    phitot = cnul;
    
    /* Electric field vector */
    for(m=0; (m<DIM); m++)
      Ei[m] = cnul;
    /* Loop over lattice vectors in fourier space */
    for(ix=0; (ix < nx); ix++) {
      for(iy=0; (iy < ny); iy++) {
	for(iz=0; (iz < nz); iz++) {
	  /* exclude k = (0,0,0) */
	  if ((ix == 0) && (iy == 0) && (iz == 0))
	    continue;	  
	  
	  calc_k(lll,ix,iy,iz,nx,ny,nz,k);

	  QQk    = qtab[ix][iy][iz];
	  Gk     = gtab[ix][iy][iz];
	  
	  /* Calculate potential */
	  phik   = phi_k(k,x[i],Gk,QQk);
	  phitot = cadd(phik,phitot);
	  
	  for(m=0; (m<DIM); m++)
	    Ei[m] = cadd(Ei[m],rcmul(k[m],phik));
	}
      }
    }
    /* Potential at atom i, and energy contribution */
#ifdef DEBUG
    fprintf(log,"phi[%3d] = %10.3e + i %10.3e\n",i,phitot.re,phitot.im);
#endif
    phi[i]  = phitot.re;
    energy += phi[i]*charge[i];
    
    /* Force on the atom */
    for(m=0; (m<DIM); m++) {
      f[i][m] = charge[i]*Ei[m].im;
#ifdef DEBUG
      if (fabs(Ei[m].re) > 1e-6)
	fprintf(log,"Ei[%1d] atom %3d = %10.3e + i%10.3e\n",
		m,i,Ei[m].re,Ei[m].im);
#endif
    }
  }
  return energy;
}
