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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_tables_c = "$Id$";

#include <math.h>
#include "typedefs.h"
#include "names.h"
#include "smalloc.h"
#include "fatal.h"
#include "futil.h"
#include "xvgr.h"
#include "vec.h"
#include "main.h"
#include "network.h"
 
/* All the possible (implemented) table functions */
enum { 
  etabLJ6,   etabLJ12, etabLJ6Shift, etabLJ12Shift, etabShift,
  etabRF,    etabCOUL, etabEwald, etabLJ6Switch, etabLJ12Switch,etabCOULSwitch, 
  etabEXPMIN,etabUSER, etabNR 
};
static char *tabnm[etabNR] = { 
  "LJ6",   "LJ12", "LJ6Shift", "LJ12Shift", "Shift",
  "RF",    "COUL", "Ewald", "LJ6Switch", "LJ12Switch","COULSwitch", 
  "EXPMIN","USER" 
};
/* This flag tells whether this is a Coulomb type funtion */
bool bCoulomb[etabNR] = { FALSE, FALSE, FALSE, FALSE, TRUE,
			  TRUE,  TRUE, TRUE, FALSE, FALSE, TRUE, 
			  FALSE, FALSE }; 

/* Index in the table that says which function to use */
enum { etiCOUL, etiLJ6, etiLJ12, etiNR };

typedef struct {
  int  nx,nx0;
  real tabscale;
  real *x,*v,*v2;
} t_tabledata;

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))

static void spline(real x[],real y[],int n,real yp1,real ypn,real y2[])
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

static void splint(real xa[],real ya[],real y2a[],
		   int n,real x,real *y,real *yp)
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

static void copy2table(int n,int offset,int stride,
		       real x[],real Vtab[],real Vtab2[],
		       real dest[],real r_zeros)
{
  int  i,nn0;
  real F,G,H,h;
    
  for(i=1; (i<n); i++) {
    h   = x[i+1]-x[i];
    F   = (Vtab[i+1]-Vtab[i]-(h*h/6.0)*(2*Vtab2[i]+Vtab2[i+1]));
    G   = (h*h/2.0)*Vtab2[i];
    H   = (h*h/6.0)*(Vtab2[i+1]-Vtab2[i]);
    nn0 = offset+i*stride;
    dest[nn0]   = Vtab[i];
    dest[nn0+1] = F;
    dest[nn0+2] = G;
    dest[nn0+3] = H;
  }
  if (r_zeros > 0.0) {
    for(i=1; (i<n); i++) {
      if (0.5*(x[i]+x[i+1]) >= r_zeros) {
	nn0 = offset+i*stride;
	dest[nn0]   = 0.0;
	dest[nn0+1] = 0.0;
	dest[nn0+2] = 0.0;
	dest[nn0+3] = 0.0;
      }
    }
  }
}

static void init_table(FILE *fp,int n,int nx0,int tabsel,
		       real tabscale,t_tabledata *td,bool bAlloc)
{
  int i;
  
  td->nx  = n;
  td->nx0 = nx0;
  td->tabscale = tabscale;
  if (bAlloc) {
    snew(td->x,td->nx+1);
    snew(td->v,td->nx+1);
    snew(td->v2,td->nx+1);
  }
  for(i=td->nx0; (i<=td->nx); i++)
    td->x[i] = i/tabscale;
}

static void read_tables(FILE *fp,char *fn,t_tabledata td[])
{
  char *libfn;
  real **yy=NULL;
  int  k,i,nx,nx0,ny,nny;
  bool bCont;
  real tabscale;

  nny = 2*etiNR+1;  
  libfn = low_libfn(fn,TRUE);
  nx  = read_xvg(libfn,&yy,&ny);
  if (ny != nny)
    fatal_error(0,"Trying to read file %s, but nr columns = %d, should be %d",
		libfn,ny,nny);
  bCont = TRUE;
  for(nx0=0; bCont && (nx0 < nx); nx0++)
    for(k=1; (k<ny); k++)
      if (yy[k][nx0] != 0)
	bCont = FALSE;
  if (nx0 == nx)
    fatal_error(0,"All elements in table %s are zero!\n",libfn);
    
  tabscale = (nx-1)/(yy[0][nx-1] - yy[0][0]);
  for(k=0; (k<etiNR); k++) {
    init_table(fp,nx,nx0,etabUSER,tabscale,&(td[k]),TRUE);
    for(i=0; (i<nx); i++) {
      td[k].x[i]  = yy[0][i];
      td[k].v[i]  = yy[2*k+1][i];
      td[k].v2[i] = yy[2*k+2][i];
    }
  }
  for(i=0; (i<ny); i++)
    sfree(yy[i]);
  sfree(yy);
  
  if (fp) 
    fprintf(fp,"Read user tables from %s with %d data points.\n"
	    "Tabscale = %g points/nm\n",libfn,nx,tabscale);
}

static void done_tabledata(t_tabledata *td)
{
  int i;
  
  if (!td)
    return;
    
  sfree(td->x);
  sfree(td->v);
  sfree(td->v2);
}

static void fill_table(t_tabledata *td,int tp,t_forcerec *fr)
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
  real expr,Vtab,Ftab,Vtab2,Ftab2;
  /* Parameters for David's function */
  real A=0,B=0,C=0,A_3=0,B_4=0;
  /* Parameters for the switching function */
  real ksw,swi,swi1,swi2,swi3;
  /* Temporary parameters */
  bool bSwitch,bShift;
  real VtabT;  
  real VtabT1;  
  real VtabT2; 
  real VtabT3;
  real ewc=fr->ewaldcoeff;
  real isp= 0.564189583547756;
   
  bSwitch = ((tp == etabLJ6Switch)    || (tp == etabLJ12Switch)    || 
	     (tp == etabCOULSwitch));
  bShift  = ((tp == etabLJ6Shift) || (tp == etabLJ12Shift) || 
	     (tp == etabShift));

  if (bCoulomb[tp]) {
    r1 = fr->rcoulomb_switch;
    rc = fr->rcoulomb;
  } 
  else {
    r1 = fr->rvdw_switch;
    rc = fr->rvdw;
  }
  if (bSwitch)
    ksw  = 1.0/pow3(rc-r1);
  else
    ksw  = 0.0;
  if (bShift) {
    if (tp == etabShift) 
      p=1;
    else if (tp == etabLJ6Shift) 
      p=6; 
    else 
      p=12;
    
    A = p * ((p+1)*r1-(p+4)*rc)/(pow(rc,p+2)*pow2(rc-r1));
    B = -p * ((p+1)*r1-(p+3)*rc)/(pow(rc,p+2)*pow3(rc-r1));
    C = 1.0/pow(rc,p)-A/3.0*pow(rc-r1,3)-B/4.0*pow4(rc-r1);
    if (tp == etabLJ6Shift) {
      A=-A;
      B=-B;
      C=-C;
    }
    A_3=A/3.0;
    B_4=B/4.0;
  }
  if (debug) { fprintf(debug,"Further\n"); fflush(debug); }
    
#ifdef DEBSW
  fp=xvgropen("switch.xvg","switch","r","s");
#endif
  for(i=td->nx0; (i<=td->nx); i++) {
    r     = td->x[i];
    r2    = r*r;
    r6    = 1.0/(r2*r2*r2);
    r12   = r6*r6;
    Vtab  = 0.0;
    Ftab  = 0.0;
    Vtab2 = 0.0;
    Ftab2 = 0.0;
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
      Vtab  = -r6;
      Ftab  = 6.0*Vtab/r;
      Vtab2 = 7.0*Ftab/r;
      Ftab2 = 8.0*Vtab2/r;
      break;
    case etabLJ6Switch:
    case etabLJ6Shift:
      /* Dispersion */
      if (r < rc) {      
	Vtab  = -r6;
	Ftab  = 6.0*Vtab/r;
	Vtab2 = 7.0*Ftab/r;
	Ftab2 = 8.0*Vtab2/r;
      }
      break;
    case etabLJ12:
      /* Repulsion */
      Vtab  = r12;
      Ftab  = 12.0*Vtab/r;
      Vtab2 = 13.0*Ftab/r;
      Ftab2 = 14.0*Vtab2/r;
      break;
    case etabLJ12Switch:
    case etabLJ12Shift:
      /* Repulsion */
      if (r < rc) {      
	Vtab  = r12;
	Ftab  = 12.0*Vtab/r;
	Vtab2 = 13.0*Ftab/r;
	Ftab2 = 14.0*Vtab2/r;
      }  
      break;
    case etabCOUL:
      Vtab  = 1.0/r;
      Ftab  = 1.0/r2;
      Vtab2 = 2.0/(r*r2);
      Ftab2 = 6.0/(r2*r2);
      break;
    case etabCOULSwitch:
    case etabShift:
      if (r < rc) { 
	Vtab  = 1.0/r;
	Ftab  = 1.0/r2;
	Vtab2 = 2.0/(r*r2);
	Ftab2 = 6.0/(r2*r2);
      }
      break;
    case etabEwald:
      Vtab  = erfc(ewc*r)/r;
      Ftab  = erfc(ewc*r)/r2+2*exp(-(ewc*ewc*r2))*ewc*isp/r;
      Vtab2 = 2*erfc(ewc*r)/(r*r2)+4*exp(-(ewc*ewc*r2))*ewc*isp/r2+
	  4*ewc*ewc*ewc*exp(-(ewc*ewc*r2))*isp;
      Ftab2 = 6*erfc(ewc*r)/(r2*r2)+
	  12*exp(-(ewc*ewc*r2))*ewc*isp/(r*r2)+
	  8*ewc*ewc*ewc*exp(-(ewc*ewc*r2))*isp/r+
	  8*ewc*ewc*ewc*ewc*ewc*r*exp(-(ewc*ewc*r2))*isp;
      break;
    case etabRF:
      Vtab  = 1.0/r      +   fr->k_rf*r2 - fr->c_rf;
      Ftab  = 1.0/r2     - 2*fr->k_rf*r;
      Vtab2 = 2.0/(r*r2) + 2*fr->k_rf;
      Ftab2 = 6.0/(r2*r2);
      break;
    case etabEXPMIN:
      expr  = exp(-r);
      Vtab  = expr;
      Ftab  = expr;
      Vtab2 = expr;
      Ftab2 = expr;
      break;
    default:
      fatal_error(0,"Table type %d not implemented yet. (%s,%d)",
		  tp,__FILE__,__LINE__);
    }
    if (bShift) {
      /* Normal coulomb with cut-off correction for potential */
      if (r < rc) {
	Vtab -= C;
	/* If in Shifting range add something to it */
	if (r > r1) {
	  r12 = (r-r1)*(r-r1);
	  r13 = (r-r1)*r12;
	  Vtab  += - A_3*r13 - B_4*r12*r12;
	  Ftab  +=   A*r12 + B*r13;
	  Vtab2 += - 2.0*A*(r-r1) - 3.0*B*r12;
	  Ftab2 +=   2.0*A + 6.0*B*(r-r1);
	}
      }
    }
    
    if ((r > r1) && bSwitch) {
      VtabT     = Vtab;
      VtabT1    = -Ftab;
      VtabT2    = Vtab2;
      VtabT3    = -Ftab2;
      Vtab   = VtabT*swi;
      Vtab2  = VtabT2*swi + VtabT1*swi1 + VtabT1*swi1 + VtabT*swi2;
      /* Ftab   = -(VtabT1*swi+ VtabT*swi1);
	 Ftab2  = -(VtabT3*swi + VtabT2*swi1 + VtabT1*swi2 + VtabT2*swi1 +
	 VtabT2*swi1 + VtabT1*swi2 + VtabT1*swi2 + VtabT*swi3);
      */
    }  
    
    /* Ftab  /= r;
       Ftab2 /= r;
    */
    td->v[i]  = Vtab;
    td->v2[i] = Vtab2;
  }

#ifdef DEBSW
  fclose(fp);
#endif
}

static void set_table_type(int tabsel[],t_forcerec *fr)
{
  /* Set the different table indices.
   * Coulomb first.
   */
  
  switch (fr->eeltype) {
  case eelCUT:
    tabsel[etiCOUL] = etabCOUL;
    break;
  case eelPPPM:
  case eelPOISSON:
    if ((fr->rcoulomb > fr->rcoulomb_switch) && fr->bcoultab)
	tabsel[etiCOUL] = etabShift;
    else
	tabsel[etiCOUL] = etabCOUL;  /* 1-4 */
    break;
  case eelSHIFT:
    if (fr->rcoulomb > fr->rcoulomb_switch)
      tabsel[etiCOUL] = etabShift;
    else
      tabsel[etiCOUL] = etabCOUL;
    break;
  case eelEWALD:
  case eelPME:
      if(fr->bcoultab)
	  tabsel[etiCOUL] = etabEwald;
      else
	  tabsel[etiCOUL] = etabCOUL; /* 1-4 */
      break;
  case eelRF:
  case eelGRF:
    tabsel[etiCOUL] = etabRF;
    break;
  case eelSWITCH:
    tabsel[etiCOUL] = etabCOULSwitch;
    break;
  case eelUSER:
    tabsel[etiCOUL] = etabUSER;
    break;
  default:
    fatal_error(0,"Invalid eeltype %d in %s line %d",fr->eeltype,
		__FILE__,__LINE__);
  }
  
  /* Van der Waals time */
  if (fr->bBHAM) {
    tabsel[etiLJ6]  = etabLJ6;
    tabsel[etiLJ12] = etabEXPMIN;
  }
  else {
    switch (fr->vdwtype) {
    case evdwSWITCH:
      tabsel[etiLJ6]  = etabLJ6Switch;
      tabsel[etiLJ12] = etabLJ12Switch;
      break;
    case evdwSHIFT:
      tabsel[etiLJ6]  = etabLJ6Shift;
      tabsel[etiLJ12] = etabLJ12Shift;
      break;
    case evdwUSER:
      tabsel[etiLJ6]  = etabUSER;
      tabsel[etiLJ12] = etabUSER;
      break;
    case evdwCUT:
      tabsel[etiLJ6]  = etabLJ6;
      tabsel[etiLJ12] = etabLJ12;
      break;
    default:
      fatal_error(0,"Invalid vdwtype %d in %s line %d",fr->vdwtype,
		  __FILE__,__LINE__);
    } 
  }
}

void make_tables(FILE *out,t_forcerec *fr,bool bVerbose,char *fn)
{
  static char *fns[3] = { "ctab.xvg", "dtab.xvg", "rtab.xvg" };
  FILE        *fp;
  t_tabledata *td;
  bool        bReadTab,bGenTab;
  real        x0,y0,yp,rtab;
  int         i,j,k,nx,nx0,tabsel[etiNR];
 
  set_table_type(tabsel,fr);
  snew(td,etiNR);
  fr->tabscale = 0;
  rtab         = fr->rtab;
  nx0          = 10;
  nx           = 0;
  
  /* Check whether we have to read or generate */
  bReadTab = FALSE;
  bGenTab  = FALSE;
  for(i=0; (i<etiNR); i++) {
    if (tabsel[i] == etabUSER) 
      bReadTab = TRUE;
    else
      bGenTab  = TRUE;
  }
  if (bReadTab) {
    read_tables(out,fn,td);
    fr->tabscale = td[0].tabscale;
    fr->rtab     = td[0].x[td[0].nx-1];
    nx0          = td[0].nx0;
    nx           = fr->ntab = fr->rtab*fr->tabscale;
    if (fr->rtab < rtab) 
      fatal_error(0,"Tables in file %s not long enough for cut-off:\n"
		  "\tshould be at least %f nm\n",fn,rtab);
  }
  if (bGenTab) {
    if (!bReadTab) {
#ifdef DOUBLE
      fr->tabscale = 2000.0;
#else
      fr->tabscale = 500.0;
#endif
      nx = fr->ntab = fr->rtab*fr->tabscale;
    }
  }
  snew(fr->coulvdwtab,12*(nx+1)+1);
  for(k=0; (k<etiNR); k++) {
    if (tabsel[k] != etabUSER) {
      init_table(out,nx,nx0,tabsel[k],
		 (tabsel[k] == etabEXPMIN) ? fr->tabscale_exp : fr->tabscale,
		 &(td[k]),!bReadTab);
      fill_table(&(td[k]),tabsel[k],fr);
      if (out) 
	fprintf(out,"Generated table with %d data points for %s.\n"
		"Tabscale = %g points/nm\n",
		td[k].nx,tabnm[tabsel[k]],td[k].tabscale);

    }
    copy2table(td[k].nx,k*4,12,td[k].x,td[k].v,td[k].v2,fr->coulvdwtab,-1);
  
    if (bDebugMode() && bVerbose) {
      fp=xvgropen(fns[k],fns[k],"r","V"); 
      for(i=td[k].nx0; (i<td[k].nx); i++) {
	for(j=0; (j<4); j++) {
	  x0=td[k].x[i]+0.25*j*(td[k].x[i+1]-td[k].x[i]);
	  splint(td[k].x,td[k].v,td[k].v2,nx-3,x0,&y0,&yp);
	  fprintf(fp,"%15.10e  %15.10e  %15.10e\n",x0,y0,yp);
	}
      }
      ffclose(fp);
    }
    done_tabledata(&(td[k]));
  }
  sfree(td);
}

