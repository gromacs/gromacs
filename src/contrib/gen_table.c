#include "math.h"
#include "stdio.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "statutil.h"

enum { mGuillot, mN61, mMaaren, mNR };

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

static void do_guillot(FILE *fp,int eel,double resolution)
{
  double qO     = -0.888;
  double qOd    = 0.226;
  double f0     = qOd/qO;
  double xi     = 0.15;
  double xir    = 0.0615;
  double sqpi   = sqrt(M_PI);
  int    i,i0,imax;
  double z,r,r1,r2,vc,vc2,vd,vd2,vr,vr2,vrep,vrep2;

  imax = 3/resolution;
  for(i=0; (i<=imax); i++) {
    r     = i*resolution;
    /* Avoid very large numbers */
    if (r < 0.04) {
      vc = vc2 = vd = vd2 = vrep = vrep2 = 0;
    }
    else {
      r1    = r/(2*xi);
      r2    = r/(sqrt(2)*xi);
      vc    = (1+sqr(f0)*erf(r1) + 2*f0*erf(r2))/r;
      vc2   = ((2/sqr(r))*(vc -
			   sqr(f0)*erf1(r1)/(2*xi) -
			   4*f0*erf1(r2)/sqrt(2)*xi) + 
	       (1/r)*(sqr(f0/(2.0*xi))*erf2(r1) + (2*f0/sqr(xi)))*erf2(r2));
      vd    = -1.0/(r*r*r*r*r*r);
      vd2   = 42.0*vd/(r*r);
      z     = r/(2.0*xir);
      vrep  = erfc(z)/z;
      vrep2 = (sqpi*vrep/(2.0*z*z)+(1.0/(z*z)+1)*exp(-z*z))/(sqpi*sqr(xir));
    }
    fprintf(fp,"%12.5e  %12.5e  %12.5e   %12.5e  %12.5e  %12.5e  %12.5e\n",
	    r,vc,vc2,vd,vd2,vrep,vrep2);
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "gen_table generates tables for mdrun for use with the USER defined",
    "potentials."
  };
  static char *opt[]     = { NULL, "cut", "rf", "pme", NULL };
  static char *model[]   = { NULL, "guillot", "n61", "maaren", NULL };
  static real resolution = 0.001;
  static int  npow       = 12;
  t_pargs pa[] = {
    { "-e",      FALSE, etENUM, {opt},
      "Electrostatics type: cut, rf or pme" },
    { "-m",      FALSE, etENUM, {model},
      "Model for the tables" },
    { "-resol",  FALSE, etREAL, {&resolution},
      "Resolution of the table (nm)" },
    { "-n",      FALSE, etINT,  {&npow},
      "Power for the repulsion potential (only with model n61)" }
  };
#define NPA asize(pa)
  t_filenm fnm[] = {
    { efXVG, "-o", "table", ffWRITE }
  };
#define NFILE asize(fnm)
  FILE *fp;
  int  eel,m;
  
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
  else 
    gmx_fatal(FARGS,"Invalid argument %s for option -m",opt[0]);
    
  fp = ffopen(opt2fn("-o",NFILE,fnm),"w");
  switch (m) {
  case mGuillot:
    do_guillot(fp,eel,resolution);
    break;
  case mN61:
    do_n61(fp,eel,resolution,npow);
    break;
  case mMaaren:
  default:
    gmx_fatal(FARGS,"Model %s not supported yet",model[0]);
  }  
  fclose(fp);
  
  thanx(stdout);
  
  return 0;
}
