/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_princ_c = "$Id$";

#include "typedefs.h"
#include "vec.h"
#include "smalloc.h"
#include "gstat.h"

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

void pr_jacobi(real **a,int n,real d[],real **v,int *nrot)
{
  int j,iq,ip,i;
  real tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
  
  snew(b,n+1);
  snew(z,n+1);
  for (ip=1;ip<=n;ip++) {
    for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=1;ip<=n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++)
	sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      sfree(z);
      sfree(b);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
	    && fabs(d[iq])+g == fabs(d[iq]))
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if (fabs(h)+g == fabs(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
	  for (j=1;j<=ip-1;j++) {
	    ROTATE(a,j,ip,j,iq)
	      }
	  for (j=ip+1;j<=iq-1;j++) {
	    ROTATE(a,ip,j,j,iq)
	      }
	  for (j=iq+1;j<=n;j++) {
	    ROTATE(a,ip,j,iq,j)
	      }
	  for (j=1;j<=n;j++) {
	    ROTATE(v,j,ip,j,iq)
	      }
	  ++(*nrot);
	}
      }
    }
    for (ip=1;ip<=n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  fatal_error(0,"Too many iterations in routine JACOBI");
}

#undef ROTATE

static void m_op(matrix mat,rvec x)
{
  rvec xp;
  int  m;
  
  for(m=0; (m<DIM); m++)
    xp[m]=mat[m][XX]*x[XX]+mat[m][YY]*x[YY]+mat[m][ZZ]*x[ZZ];
  fprintf(stderr,"x    %8.3f  %8.3f  %8.3f\n",x[XX],x[YY],x[ZZ]);
  fprintf(stderr,"xp   %8.3f  %8.3f  %8.3f\n",xp[XX],xp[YY],xp[ZZ]);
  fprintf(stderr,"fac  %8.3f  %8.3f  %8.3f\n",xp[XX]/x[XX],xp[YY]/x[YY],
	  xp[ZZ]/x[ZZ]);
}

#define NDIM 4

static void ptrans(char *s,real **inten,real d[],real e[])
{
  int  m;
  real n,x,y,z;
#ifdef DEBUG  
  for(m=1; (m<NDIM); m++) {
    x=inten[m][1];
    y=inten[m][2];
    z=inten[m][3];
    n=x*x+y*y+z*z;
    fprintf(stderr,"%8s %8.3f %8.3f %8.3f, norm:%8.3f, d:%8.3f, e:%8.3f\n",
	    s,x,y,z,sqrt(n),d[m],e[m]);
  }
  fprintf(stderr,"\n");
#endif
}

void t_trans(matrix trans,real d[],real **ev)
{
  rvec x;
  int  j;
#ifdef DEBUG  
  for(j=0; (j<DIM); j++) {
    x[XX]=ev[1][j+1];
    x[YY]=ev[2][j+1];
    x[ZZ]=ev[3][j+1];
    m_op(trans,x);
    fprintf(stderr,"d[%d]=%g\n",j,d[j+1]);
  }
#endif
}

void principal_comp(int n,atom_id index[],t_atom atom[],rvec x[],
		    matrix trans,rvec d)
{
  int  i,j,ai,m,nrot;
  real mm,rx,ry,rz;
  real **inten,dd[NDIM],e[NDIM],tvec[NDIM],**ev;
  real temp;
  
  snew(inten,NDIM);
  snew(ev,NDIM);
  for(i=0; (i<NDIM); i++) {
    snew(inten[i],NDIM);
    snew(ev[i],NDIM);
    dd[i]=e[i]=0.0;
  }
    
  for(i=0; (i<NDIM); i++)
    for(m=0; (m<NDIM); m++)
      inten[i][m]=0;
  for(i=0; (i<n); i++) {
    ai=index[i];
    mm=atom[ai].m;
    rx=x[ai][XX];
    ry=x[ai][YY];
    rz=x[ai][ZZ];
    inten[1][1]+=mm*(sqr(ry)+sqr(rz));
    inten[2][2]+=mm*(sqr(rx)+sqr(rz));
    inten[3][3]+=mm*(sqr(rx)+sqr(ry));
    inten[2][1]-=mm*(ry*rx);
    inten[3][1]-=mm*(rx*rz);
    inten[3][2]-=mm*(rz*ry);
  }
  inten[1][2]=inten[2][1];
  inten[1][3]=inten[3][1];
  inten[2][3]=inten[3][2];
  ptrans("initial",inten,dd,e);
  
  for(i=0; (i<DIM); i++) {
    for(m=0; (m<DIM); m++)
      trans[i][m]=inten[1+i][1+m];
  }

  /* Call numerical recipe routines */
  pr_jacobi(inten,3,dd,ev,&nrot);
  ptrans("jacobi",ev,dd,e);
  
  /* Sort eigenvalues in descending order */
#define SWAPPER(i) 			\
  if (fabs(dd[i+1]) > fabs(dd[i])) {	\
    temp=dd[i];			\
    for(j=0; (j<NDIM); j++) tvec[j]=ev[j][i];\
    dd[i]=dd[i+1];			\
    for(j=0; (j<NDIM); j++) ev[j][i]=ev[j][i+1];		\
    dd[i+1]=temp;			\
    for(j=0; (j<NDIM); j++) ev[j][i+1]=tvec[j];			\
  }
  SWAPPER(1);
  SWAPPER(2);
  SWAPPER(1);
  ptrans("swap",ev,dd,e);
  
  t_trans(trans,dd,ev);
    
  for(i=0; (i<DIM); i++) {
    d[i]=dd[i+1];
    for(m=0; (m<DIM); m++)
      trans[i][m]=ev[1+m][1+i];
  }
    
  for(i=0; (i<NDIM); i++) {
    sfree(inten[i]);
    sfree(ev[i]);
  }
  sfree(inten);
  sfree(ev);
}

void rotate_atoms(int gnx,atom_id index[],rvec x[],matrix trans)
{
  real   xt,yt,zt;
  int    i,ii;
  
  for(i=0; (i<gnx); i++) {
    ii=index[i];
    xt=x[ii][XX];
    yt=x[ii][YY];
    zt=x[ii][ZZ];
    x[ii][XX]=trans[XX][XX]*xt+trans[XX][YY]*yt+trans[XX][ZZ]*zt;
    x[ii][YY]=trans[YY][XX]*xt+trans[YY][YY]*yt+trans[YY][ZZ]*zt;
    x[ii][ZZ]=trans[ZZ][XX]*xt+trans[ZZ][YY]*yt+trans[ZZ][ZZ]*zt;
  }
}

real calc_xcm(rvec x[],int gnx,atom_id index[],t_atom atom[],rvec xcm,
	      bool bQ)
{
  int  i,ii,m;
  real m0,tm;

  clear_rvec(xcm);
  tm=0;
  for(i=0; (i<gnx); i++) {
    ii=index[i];
    if (bQ)
      m0=fabs(atom[ii].q);
    else
      m0=atom[ii].m;
    tm+=m0;
    for(m=0; (m<DIM); m++)
      xcm[m]+=m0*x[ii][m];
  }
  for(m=0; (m<DIM); m++)
    xcm[m]/=tm;
  
  return tm;
}

real sub_xcm(rvec x[],int gnx,atom_id index[],t_atom atom[],rvec xcm,
	     bool bQ)
{
  int  i,ii;
  real tm;
  
  tm=calc_xcm(x,gnx,index,atom,xcm,bQ);
  for(i=0; (i<gnx); i++) {
    ii=index[i];
    rvec_dec(x[ii],xcm);
  }
  return tm;
}

void add_xcm(rvec x[],int gnx,atom_id index[],rvec xcm)
{
  int  i,ii;
  
  for(i=0; (i<gnx); i++) {
    ii=index[i];
    rvec_inc(x[ii],xcm);
  }
}

