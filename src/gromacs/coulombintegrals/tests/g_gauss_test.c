#include <stdio.h>
#include "gaussian_integrals.h"

int main(int argc,char *argv[])
{
  FILE *fp;
  int  i,n=1000;
  double xi;
  double r;
  
  fp = fopen("test_gauss.xvg","w");
  i = 0;
  for(xi=20; (xi<=100); xi+= 20) 
    fprintf(fp,"@ s%d legend \"Nuclear xi = %g\"\n",i++,xi);
  for(xi=20; (xi<=100); xi+= 20) 
    fprintf(fp,"@ s%d legend \"Coulomb xi = %g  xj = %g\"\n",i++,xi,xi);
  for(xi=20; (xi<=100); xi+= 20) {
    fprintf(fp,"@type xy\n");
    for(i=0; (i<n); i++) {
      r = i*1.0/n;
      fprintf(fp,"%12e  %12e\n",r,Nuclear_GG(r,xi));
    }
    fprintf(fp,"&\n");
  }
  for(xi=20; (xi<=100); xi+= 20) {
    fprintf(fp,"@type xy\n");
    for(i=0; (i<n); i++) {
      r = i*1.0/n;
      fprintf(fp,"%12e  %12e\n",r,Coulomb_GG(r,xi,xi));
    }
    fprintf(fp,"&\n");
  }
  fclose(fp);
  
  return 0;
}
