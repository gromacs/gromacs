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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_tables_c = "$Id$";

#include <math.h>
#include "typedefs.h"
#include "smalloc.h"
#include "fatal.h"
#include "futil.h"
#include "xvgr.h"
#include "vec.h"

void spline(real x[],real y[],int n,real yp1,real ypn,real y2[])
{
  int  i,k;
  real p,qn,sig,un,*u;
  
  snew(u,n+1);
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  sfree(u);
}

void splint(real xa[],real ya[],real y2a[],int n,real x,real *y,real *yp)
{
  int  klo,khi,k;
  real h,b,a,eps;
  
  klo=1;
  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) fatal_error(0,"Bad XA input to routine SPLINT");
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  *yp=(ya[khi]-ya[klo])/h+((3*a*a-1)*y2a[klo]-(3*b*b-1)*y2a[khi])*h/6.0;
  
  {
    real F,G,H;
    
    eps = b;
    F = (ya[khi]-ya[klo]-(h*h/6.0)*(2*y2a[klo]+y2a[khi]));
    G = (h*h/2.0)*y2a[klo];
    H   = (h*h/6.0)*(y2a[khi]-y2a[klo]);
    *y  = ya[klo] + eps*F + eps*eps*G + eps*eps*eps*H;
    *yp = (F + 2*eps*G + 3*eps*eps*H)/h;
  }
}

void copy2table(int n,int n0,int stride,
		real x[],real Vtab[],real Vtab2[],real dest[],real r_zeros)
{
  int  i;
  real F,G,H,h;
    
  for(i=1; (i<n); i++) {
    h = x[i+1]-x[i];
    F = (Vtab[i+1]-Vtab[i]-(h*h/6.0)*(2*Vtab2[i]+Vtab2[i+1]));
    G = (h*h/2.0)*Vtab2[i];
    H = (h*h/6.0)*(Vtab2[i+1]-Vtab2[i]);
    dest[n0+i*stride]   = Vtab[i];
    dest[n0+i*stride+1] = F;
    dest[n0+i*stride+2] = G;
    dest[n0+i*stride+3] = H;
  }
  if (r_zeros > 0.0) {
    for(i=1; (i<n); i++) {
      if (0.5*(x[i]+x[i+1]) >= r_zeros) {
	dest[n0+i*stride]   = 0.0;
	dest[n0+i*stride+1] = 0.0;
	dest[n0+i*stride+2] = 0.0;
	dest[n0+i*stride+3] = 0.0;
      }
    }
  }
}

enum { etabLJ6, etabLJ12, etabLJ6David, etabLJ12David, etabDavid,
       etabRF, etabCOUL, etabLJ6sw, etabLJ12sw, etabCOULsw, etabNR };

void read_table(int n0,int n,real x[],
		real Vtab[],real Vtab2[],
		real Ftab[],real Ftab2[],
		char *fn,t_forcerec *fr)
{
  real **y=NULL;
  int  i,nx,ny;
  
  nx = read_xvg(fn,&y,&ny);
  if (ny != 5)
    fatal_error(0,"Trying to read file %s, but no colums = %d, should be 5",
		fn,ny);
  for(i=0; (i<nx); i++) {
    x[i]     = y[0][i];
    Vtab[i]  = y[1][i];
    Ftab[i]  = y[2][i];
    Vtab2[i] = y[3][i];
    Ftab2[i] = y[4][i];
  }
  for(i=0; (i<ny); i++)
    sfree(y[i]);
  sfree(y);
}

void fill_table(int n0,int n,real x[],
		real Vtab[],real Vtab2[],
		real Ftab[],real Ftab2[],
		int tp,t_forcerec *fr)
{
  /* Calculate potential and 2nd derivative and Force and
   * second derivative!
   */
#ifdef DEBSW
  FILE *fp;
#endif
  int  i,p;
  real r1,rc,r12,r13;
  real r,r2,r6;
  real k_rf,c_rf,rffac2;
  /* Parameters for David's function */
  real A,B,C,A_3,B_4;
  /* Parameters for the switching function */
  real ksw,swi,swi1,swi2,swi3;
  /* Temporary parameters */
  bool bSwitch,bDavid;
  real VtabT;  
  real VtabT1;  
  real VtabT2; 
  real VtabT3;
   
  bSwitch= ((tp == etabLJ6sw) || (tp == etabLJ12sw) || (tp == etabCOULsw));
  bDavid= ((tp == etabLJ6David) || (tp == etabLJ12David) || (tp == etabDavid));
  r1     = fr->r1;
  rc     = fr->rc;
  k_rf   = fr->k_rf;
  c_rf   = fr->c_rf;
  rffac2 = k_rf*2.0;
  if (bSwitch)
    ksw  = 1.0/pow((rc-r1),3.0);
  else
    ksw  = 0.0;
  if (bDavid) {
    if (tp==etabDavid) 
      p=1;
    else
      if (tp==etabLJ6David) 
	p=6; 
      else 
	p=12;
    A = p * ((p+1)*r1-(p+4)*rc)/(pow(rc,p+2)*pow(rc-r1,2));
    B = -p * ((p+1)*r1-(p+3)*rc)/(pow(rc,p+2)*pow(rc-r1,3));
    C = pow(rc,-p)-A/3.0*pow(rc-r1,3)-B/4.0*pow(rc-r1,4);
    if (tp==etabLJ6David) {
      A=-A;
      B=-B;
      C=-C;
    }
    A_3=A/3.0;
    B_4=B/4.0;
  }
    
#ifdef DEBSW
  fp=xvgropen("switch.xvg","switch","r","s");
#endif
  for(i=n0; (i<=n); i++) {
    r        = x[i];
    r2       = r*r;
    r6       = pow(r2,-3.0);
    r12      = pow(r2,-6.0);
    Vtab[i]  = 0.0;
    Ftab[i]  = 0.0;
    Vtab2[i] = 0.0;
    Ftab2[i] = 0.0;
    if (bSwitch) {
      swi      = (rc-r)*(rc-r)*(rc+2*r-3*r1)*ksw;
      swi1     = 6*(rc-r)*(r1-r)*ksw;
      swi2     = -6*(r1+rc-2*r)*ksw;
      swi3     = 12*ksw;
    }
    else {
      swi = swi1 = swi2 = swi3 = 1.0;
    }
#ifdef DEBSW
    fprintf(fp,"%10g  %10g  %10g  %10g  %10g\n",r,swi,swi1,swi2,swi3);
#endif
    switch (tp) {
    case etabLJ6:
      /* Dispersion */
      Vtab[i]  = -r6;
      Ftab[i]  = 6.0*Vtab[i]/r;
      Vtab2[i] = 7.0*Ftab[i]/r;
      Ftab2[i] = 8.0*Vtab2[i]/r;
      break;
    case etabLJ6sw:
    case etabLJ6David:
      /* Dispersion */
      if (r < rc) {      
	Vtab[i]  = -r6;
	Ftab[i]  = 6.0*Vtab[i]/r;
	Vtab2[i] = 7.0*Ftab[i]/r;
	Ftab2[i] = 8.0*Vtab2[i]/r;
      }
      break;
    case etabLJ12:
      /* Repulsion */
      Vtab[i]  = r12;
      Ftab[i]  = 12.0*Vtab[i]/r;
      Vtab2[i] = 13.0*Ftab[i]/r;
      Ftab2[i] = 14.0*Vtab2[i]/r;
      break;
    case etabLJ12sw:
    case etabLJ12David:
      /* Repulsion */
      if (r < rc) {      
	Vtab[i]  = r12;
	Ftab[i]  = 12.0*Vtab[i]/r;
	Vtab2[i] = 13.0*Ftab[i]/r;
	Ftab2[i] = 14.0*Vtab2[i]/r;
      }  
      break;
    case etabCOUL:
      Vtab[i]  = 1.0/r;
      Ftab[i]  = 1.0/r2;
      Vtab2[i] = 2.0/(r*r2);
      Ftab2[i] = 6.0/(r2*r2);
      break;
    case etabCOULsw:
    case etabDavid:
      if (r < rc) { 
	Vtab[i]  = 1.0/r;
	Ftab[i]  = 1.0/r2;
	Vtab2[i] = 2.0/(r*r2);
	Ftab2[i] = 6.0/(r2*r2);
      }
      break;
    case etabRF:
      Vtab[i]  = 1.0/r      + k_rf*r2  - c_rf;
      Ftab[i]  = 1.0/r2     - rffac2*r;
      Vtab2[i] = 2.0/(r2*r) + rffac2;
      Ftab2[i] = 6.0/(r2*r2);
      break;
    default:
      fatal_error(0,"Table type %d not implemented yet. (%s,%d)",
		  tp,__FILE__,__LINE__);
    }
    if (bDavid) {
      /* Normal coulomb with cut-off correction for potential */
      if (r < rc) {
	Vtab[i] -= C;
	/* If in shifting range add something to it */
	if (r > r1) {
	  r12 = (r-r1)*(r-r1);
	  r13 = (r-r1)*r12;
	  Vtab[i]  += - A_3*r13 - B_4*r12*r12;
	  Ftab[i]  +=   A*r12 + B*r13;
	  Vtab2[i] += - 2.0*A*(r-r1) - 3.0*B*r12;
	  Ftab2[i] +=   2.0*A + 6.0*B*(r-r1);
	}
      }
    }
    
    if ((r > r1) && bSwitch) {
      VtabT     = Vtab[i];
      VtabT1    = -Ftab[i];
      VtabT2    = Vtab2[i];
      VtabT3    = -Ftab2[i];
      Vtab[i]   = VtabT*swi;
      Ftab[i]   = -(VtabT1*swi+ VtabT*swi1);
      Vtab2[i]  = VtabT2*swi + VtabT1*swi1 + VtabT1*swi1 + VtabT*swi2;
      Ftab2[i]  = -(VtabT3*swi + VtabT2*swi1 + VtabT1*swi2 + VtabT2*swi1 +
		    VtabT2*swi1 + VtabT1*swi2 + VtabT1*swi2 + VtabT*swi3);
    }  

    Ftab[i]  /= r;
    Ftab2[i] /= r;
  }
#ifdef DEBSW
  fclose(fp);
#endif
}

void make_tables(t_forcerec *fr,bool bVerbose)
{
  char *fns[3]    = { "ctab.xvg", "dtab.xvg", "rtab.xvg"  };
  /*
  int  pppmtab[3] = { etabDavid,  etabLJ6David,   etabLJ12David };
  int  rftab[3]   = { etabRF,     etabLJ6,   etabLJ12 };
  int  normtab[3] = { etabCOUL,   etabLJ6,   etabLJ12 };
  int  swtab[3]   = { etabCOULsw, etabLJ6sw, etabLJ12sw  };
  */
  int  ljtab[2]   = { etabLJ6,      etabLJ12 };
  int  ljshtab[2] = { etabLJ6David, etabLJ12David };
  int  ljswtab[2] = { etabLJ6sw,    etabLJ12sw  };
  FILE     *fp;
  real     x0,y0,yp;
  int      i,j,k,n0,n,tabsel;
  real     *x,*Vtab,*Vtab2,*Ftab,*Ftab2;
  
#ifdef DOUBLE
  fr->tabscale = 2000.0;
#else
  fr->tabscale = 500.0;
#endif
  n = fr->ntab = (fr->rc+0.5)*fr->tabscale;

  snew(fr->VFtab,12*n+1);
  snew(x,n+1);
  
  n0 = 10;
  for(i=n0; (i<=n); i++) {
    x[i]  = (i/fr->tabscale);
  }
  
  /* Now fill three tables with 
   * Dispersion, Repulsion and Coulomb (David's function)
   * or Reaction Field, or Plain Coulomb
   * respectively.
   */
  /* Temp storage */
  snew(Vtab,n+1);
  snew(Vtab2,n+1);
  snew(Ftab,n+1);
  snew(Ftab2,n+1);

  for(k=0; (k<3); k++) {
    /* Check which table we have to use */
    if (k==0)
      switch (fr->eeltype) {
      case eelPPPM:
      case eelPOISSON:
      case eelSHIFT:
	tabsel = etabDavid;
	break;
      case eelRF:
      case eelGRF:
	tabsel = etabRF;
	break;
      case eelSWITCH:
	tabsel = etabCOULsw;
	break;
      case eelUSER:
      default:
	tabsel = etabCOUL;
      }
    else
      switch (fr->eeltype) {
      case eelSWITCH:
	tabsel = ljswtab[k-1];
	break;
      default:
	if (fr->bLJshift)
	  tabsel = ljshtab[k-1];
	else
	  tabsel = ljtab[k-1];
      }  
    if (fr->eeltype == eelUSER) {
      /* Read tables from file */
      read_table(n0,n,x,Vtab,Vtab2,Ftab,Ftab2,fns[k],fr);
    }
    else {
      fill_table(n0,n,x,Vtab,Vtab2,Ftab,Ftab2,tabsel,fr);
    }
    /*
       if (tabsel == etabRF)
       copy2table(n,k*4,12,x,Vtab,Vtab2,fr->VFtab,fr->rc);
       else
       */
    copy2table(n,k*4,12,x,Vtab,Vtab2,fr->VFtab,-1);
    
    if (bDebugMode()) {
      fp=xvgropen(fns[k],fns[k],"r","V");  
      for(i=n0; (i<n); i++) {
	for(j=0; (j<4); j++) {
	  x0=x[i]+0.25*j*(x[i+1]-x[i]);
	  splint(x,Vtab,Vtab2,n-3,x0,&y0,&yp);
	  fprintf(fp,"%15.10e  %15.10e  %15.10e\n",x0,y0,yp);
	}
      }
      ffclose(fp);
    }
  }
  sfree(Vtab);
  sfree(Vtab2);
  sfree(Ftab);
  sfree(Ftab2);
  sfree(x);
}






