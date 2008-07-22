#include <iostream>
#include <fstream>
using namespace std;
#include "slater_S_integrals.h"

#define CAL2JOULE	 (4.184)		/* id		*/
#define FACEL		 (332.0636*CAL2JOULE)   /* (10 * (ONE_4PI_EPS0)) */
#define ONE_4PI_EPS0	 (FACEL*0.1)            /* 1/(4*pi*e0)*/

static void xvg_header(FILE *fp,int bSS,double xii,double xij)
{
  int j,k,n;
  
  fprintf(fp,"@g0 on\n");
  fprintf(fp,"@with g0\n");
  fprintf(fp,"@ world 0, 0, 1, 600.0\n");
  fprintf(fp,"@    stack world 0, 0, 0, 0\n");
  fprintf(fp,"@    znorm 1\n");

  if (bSS) 
    fprintf(fp,"@ subtitle \"\\f{12}x\\f{4}\\si\\N = %g/nm  \\f{12}x\\f{4}\\sj\\N = %g/nm\"\n",xii,xij);
  else
    fprintf(fp,"@ subtitle \"\\f{12}x\\f{4}\\si\\N = %g/nm\"\n",xii);
  fprintf(fp,"@    subtitle font 4\n");
  fprintf(fp,"@    subtitle size 1.500000\n");

  fprintf(fp,"@ xaxis label \"r (nm)\"\n");
  fprintf(fp,"@ xaxis  label font 4\n");
  fprintf(fp,"@ xaxis  label char size 1.500000\n");
  fprintf(fp,"@ xaxis  ticklabel char size 1.250000\n");
  fprintf(fp,"@ xaxis  ticklabel font 4\n");
  fprintf(fp,"@    xaxis  tick major 0.2\n");
  fprintf(fp,"@    xaxis  tick minor ticks 1\n");

  fprintf(fp,"@ yaxis label \"E\\sCoulomb\\N kJ/mol e\\S2\"\n"); 
  fprintf(fp,"@ yaxis  label font 4\n");
  fprintf(fp,"@ yaxis  label char size 1.500000\n");
  fprintf(fp,"@ yaxis  ticklabel char size 1.250000\n");
  fprintf(fp,"@ yaxis  ticklabel font 4\n");
  fprintf(fp,"@    yaxis  tick major 100.0\n");
  fprintf(fp,"@    yaxis  tick minor ticks 3\n");
  fprintf(fp,"@    legend font 4\n");
  fprintf(fp,"@    legend char size 1.500000\n");

  n = 0;
  if (0) {
    fprintf(fp,"@ s%d legend \"1/r\"\n",n);
    fprintf(fp,"@ s%d line linewidth 2.0\n",n);
    n++;
  }
  if (bSS) {
    for(j=1; (j<=SLATER_MAX); j++) {
      for(k=j; (k<=SLATER_MAX); k++) {
	fprintf(fp,"@ s%d legend \"%dS - %dS\"\n",n,k,j);
	fprintf(fp,"@ s%d line linewidth 2.0\n",n);
	n++;
      }
    }
  }
  else {
    for(j=1; (j<=SLATER_MAX); j++) {
      fprintf(fp,"@ s%d legend \"%dS\"\n",n,j);
      fprintf(fp,"@ s%d line linewidth 2.0\n",n);
      n++;
    }
  }
}

static void dump_ss(char *fn,cl_F xi,cl_F xj,int nmax)
{
  int      i,j,k;
  char     buf[256];
  cl_F     r,dr,S;
  double   rd,Sd;
  ofstream fp;
  
  sprintf(buf,"0.0_%d",PRECISION);
  S = buf;
  sprintf(buf,"0.001_%d",PRECISION);
  r = dr = buf;
  if (0) {
    fp.open(fn,fstream::app);
    fp << "@type xy\n";
    for(i=1; (i<=nmax); i++) {
      S = ONE_4PI_EPS0/r;
      rd = double_approx(r);
      Sd = double_approx(S);
      fp << scientific << rd << "  " << Sd << endl;
      fp.flush();
      r += dr;
    }
    fp << "&" << endl;
    fp.close();
  }
  for(j=1; (j<=SLATER_MAX); j++) {
    for(k=j; (k<=SLATER_MAX); k++) {
      cout << "Dumping " << j << " " << k << endl;
      fp.open(fn,fstream::app);
      fp << "@type xy" << endl;
      r  = dr;
      for(i=0; (i<=nmax); i++) {
	if ((j == 1) && (k == 1))
	  S  = Slater_1S_1S(r,xi,xj);
	else
	  S  = Slater_SS[j-1][k-1](r,xi,xj);
	rd = double_approx(r);
	Sd = double_approx(S);
	fp << rd << "  " << Sd << endl;
	fp.flush();
	r += dr;
      }
      fp << "&" << endl;
      fp.close();
    }
  }
}

static void dump_ns(char *fn,cl_F xi,int nmax)
{
  int     i,j;
  cl_F r,S;
  ofstream fp;
  
  fp.open(fn,fstream::app);
  if (0) {
    fp << "type xy\n";
    for(i=1; (i<=nmax); i++) {
      r = (i*1.0/nmax);
      S = ONE_4PI_EPS0/r;
      fp << scientific << r << "  " << S << endl;
    }
    fp << "&\n";
  }
  for(j=1; (j<=SLATER_MAX); j++) {
    fp << "@type xy\n";
    for(i=1; (i<=nmax); i++) {
      r = (i*1.0/nmax);
      S = Slater_NS[j-1](r,xi);
      fp << scientific << r << "  " << S*ONE_4PI_EPS0 << endl;
    }
    fp << "&\n";
  }
  fp.close();
}
  
int main(int argc,char *argv[])
{
  FILE    *fp;
  char    buf[256],fn[256];
  cl_F    xi,xj;
  double  dxi,dxj;
  int     nmax=1000;
  
  if (argc < 3) { 
    cerr << "Usage: coulomb xi xj\n" << endl;
    exit(1);
  }
  sprintf(buf,"%s_%d",argv[1],PRECISION);
  xi = buf;
  sprintf(buf,"%s_%d",argv[2],PRECISION);
  xj = buf;
  dxi = atof(argv[1]);
  dxj = atof(argv[2]);
    
  sprintf(fn,"slater_%g_%g.xvg",dxi,dxj);
  fp = fopen(fn,"w");
  xvg_header(fp,1,dxi,dxj);
  fclose(fp);
  dump_ss(fn,xi,xj,nmax);
  
  sprintf(fn,"nuclear_%g.xvg",dxi);
  fp = fopen(fn,"w");
  xvg_header(fp,0,dxi,dxj);
  fclose(fp);
  dump_ns(fn,xi,nmax);
  
  return 0;
}
