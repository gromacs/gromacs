#include "math.h"
#include "string.h"
#include "copyrite.h"
#include "stdio.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "statutil.h"

enum { mGuillot, mN61, mMaaren, mGuillot_maple, mNR };

void do_n61(FILE *fp,int eel,double resolution,int npow)
{
  int    i,k,imax;
  double myfac[3] = { 1, -1, 1 };
  double myexp[3] = { 1, 6, 0 };
  double x,v,v2;
  
  myexp[2] = npow;
  imax     = 3.0/resolution;
  for(i=0; (i<=imax); i++) {
    x   =  i*resolution;
    
    fprintf(fp,"%10g",x);
    
    for(k=0; (k<3); k++) {
      if (x < 0.04) {
	/* Avoid very high numbers */
	v = v2 = 0;
      }
      else {
	v  =  myfac[k]*pow(x,-myexp[k]);
	v2 = (myexp[k]+1)*(myexp[k])*v/(x*x); 
      }
      fprintf(fp,"  %10g  %10g",v,v2);
    }
    fprintf(fp,"\n");
  }
}

static double erf2(double x)
{
  return -(4*x/(sqrt(M_PI)))*exp(-x*x);
}

static double erf1(double x)
{
  return (2/sqrt(M_PI))*exp(-x*x);
}

void lo_do_guillot(double r,double xi,
			    double *vc,double *vc2,
			    double *vd,double *vd2,
			    double *vr,double *vr2)
{
  double qO     = -0.888;
  double qOd    = 0.226;
  double f0     = qOd/qO;
  double xir    = 0.0615;
  double sqpi   = sqrt(M_PI);
  double r1,r2,z;
  
  r1    = r/(2*xi);
  r2    = r/(sqrt(2)*xi);
  *vc   = (1+sqr(f0)*erf(r1) + 2*f0*erf(r2))/r;
  *vc2  = ((2/sqr(r))*(*vc -
		       sqr(f0)*erf1(r1)/(2*xi) -
		       4*f0*erf1(r2)/sqrt(2)*xi) + 
	   (1/r)*(sqr(f0/(2.0*xi))*erf2(r1) + (2*f0/sqr(xi)))*erf2(r2));
  *vd   = -1.0/(r*r*r*r*r*r);
  *vd2  = 42.0*(*vd)/(r*r);
  z     = r/(2.0*xir);
  *vr   = erfc(z)/z;
  *vr2  = (sqpi*(*vr)/(2.0*z*z)+(1.0/(z*z)+1)*exp(-z*z))/(sqpi*sqr(xir));
}

void lo_do_guillot_maple(double r,double xi,
				  double *vc,double *vc2,
				  double *vd,double *vd2,
				  double *vr,double *vr2)
{
  double qO     = -0.888;
  double qOd    = 0.226;
  double f0     = qOd/qO;
  double xir    = 0.0615;
  double sqpi   = sqrt(M_PI);
  double r1,r2,z;
  
  r1    = r/(2*xi);
  r2    = r/(sqrt(2)*xi);
  *vc   = 1/r + f0*f0*erf(r/xi/2.0)/r + 2.0*f0*erf(r*sqrt(2.0)/xi/2.0)/r;
  *vc2  = 2.0/(r*r*r)-f0*f0/sqrt(M_PI)/(xi*xi*xi)*exp(-r*r/(xi
*xi)/4.0)/2.0-2.0*f0*f0/sqrt(M_PI)*exp(-r*r/(xi*xi)/4.0)/xi/(r*
r)+2.0*f0*f0*erf(r/xi/2.0)/(r*r*r)-2.0*f0/sqrt(M_PI)/(xi*xi*xi)
*exp(-r*r/(xi*xi)/2.0)*sqrt(2.0)-4.0*f0/sqrt(M_PI)*exp(-r*r/(xi
*xi)/2.0)*sqrt(2.0)/xi/(r*r)+4.0*f0*erf(r*sqrt(2.0)/xi/2.0)/(r*r*r);
  *vd   = -1.0/(r*r*r*r*r*r);
  *vd2  = -42.0/(r*r*r*r*r*r*r*r);
  z     = r/(2.0*xir);
  *vr   = 2.0*erfc(r/xir/2.0)/r*xir;
  *vr2  = 1.0/sqrt(M_PI)/(xir*xir)*exp(-r*r/(xir*xir)/4.0)+4.0/sqrt(M_PI)*exp(-r*r/(xir*xir)/4.0)/(r*r)+4.0*erfc(r/xir/2.0)/(r*r*r)*xir;
}

static void do_guillot(FILE *fp,int eel,double resolution)
{
  int    i,i0,imax;
  double xi     = 0.15;
  double r,vc,vc2,vd,vd2,vr,vr2;

  imax = 3/resolution;
  for(i=0; (i<=imax); i++) {
    r     = i*resolution;
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = vc2 = vd = vd2 = vr = vr2 = 0;
    }
    else 
      lo_do_guillot(r,xi,&vc,&vc2,&vd,&vd2,&vr,&vr2);
    fprintf(fp,"%12.5e  %12.5e  %12.5e   %12.5e  %12.5e  %12.5e  %12.5e\n",
	    r,vc,vc2,vd,vd2,vr,vr2);
  }
}

static void do_guillot_maple(FILE *fp,int eel,double resolution)
{
  int    i,i0,imax;
  double xi     = 0.15;
  double r,vc,vc2,vd,vd2,vr,vr2;

  imax = 3/resolution;
  for(i=0; (i<=imax); i++) {
    r     = i*resolution;
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = vc2 = vd = vd2 = vr = vr2 = 0;
    }
    else 
      lo_do_guillot_maple(r,xi,&vc,&vc2,&vd,&vd2,&vr,&vr2);
    fprintf(fp,"%12.5e  %12.5e  %12.5e   %12.5e  %12.5e  %12.5e  %12.5e\n",
	    r,vc,vc2,vd,vd2,vr,vr2);
  }
} 

static void do_maaren(FILE *fp,int eel,double resolution,int npow)
{
  int    i,i0,imax;
  double xi     = 0.05;
  double r,vc,vc2,vd,vd2,vr,vr2;

  imax = 3/resolution;
  for(i=0; (i<=imax); i++) {
    r     = i*resolution;
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = vc2 = vd = vd2 = vr = vr2 = 0;
    }
    else {
      lo_do_guillot(r,xi,&vc,&vc2,&vd,&vd2,&vr,&vr2);
      vr  =  pow(r,-1.0*npow);
      vr2 = (npow+1.0)*(npow)*vr/sqr(r); 
    }
    fprintf(fp,"%12.5e  %12.5e  %12.5e   %12.5e  %12.5e  %12.5e  %12.5e\n",
	    r,vc,vc2,vd,vd2,vr,vr2);
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "gen_table generates tables for mdrun for use with the USER defined",
    "potentials."
  };
  static char *opt[]     = { NULL, "cut", "rf", "pme", NULL };
  static char *model[]   = { NULL, "guillot", "n61", "maaren", "guillot_maple", NULL };
  static real resolution = 0.001;
  static int  npow       = 12;
  t_pargs pa[] = {
    { "-el",      FALSE, etENUM, {opt},
      "Electrostatics type: cut, rf or pme" },
    { "-m",      FALSE, etENUM, {model},
      "Model for the tables" },
    { "-resol",  FALSE, etREAL, {&resolution},
      "Resolution of the table (nm)" },
    { "-n",      FALSE, etINT,  {&npow},
      "Power for the repulsion potential (with model n61 or maaren)" }
  };
#define NPA asize(pa)
  t_filenm fnm[] = {
    { efXVG, "-o", "table", ffWRITE }
  };
#define NFILE asize(fnm)
  FILE *fp;
  int  eel=0,m=0;
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);
  
  if (strcmp(opt[0],"cut") == 0) 
    eel = eelCUT;
  else if (strcmp(opt[0],"rf") == 0) 
    eel = eelRF;
  else if (strcmp(opt[0],"pme") == 0) 
    eel = eelPME;
  else 
    gmx_fatal(FARGS,"Invalid argument %s for option -e",opt[0]);
  if (strcmp(model[0],"maaren") == 0) 
    m = mMaaren;
  else if (strcmp(model[0],"n61") == 0) 
    m = mN61;
  else if (strcmp(model[0],"guillot") == 0) 
    m = mGuillot;
  else if (strcmp(model[0],"guillot_maple") == 0) 
    m = mGuillot_maple;
  else 
    gmx_fatal(FARGS,"Invalid argument %s for option -m",opt[0]);
    
  fp = ffopen(opt2fn("-o",NFILE,fnm),"w");
  switch (m) {
  case mGuillot:
    do_guillot(fp,eel,resolution);
    break;
  case mGuillot_maple:
    do_guillot_maple(fp,eel,resolution);
    break;
  case mMaaren:
    do_maaren(fp,eel,resolution,npow);
    break;
  case mN61:
    do_n61(fp,eel,resolution,npow);
    break;
  default:
    gmx_fatal(FARGS,"Model %s not supported yet",model[0]);
  }  
  fclose(fp);
  
  thanx(stdout);
  
  return 0;
}
