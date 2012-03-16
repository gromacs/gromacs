#include <stdio.h>
#include <stdlib.h>
#include "slater_integrals.h"

int main(int argc,char *argv[])
{
  double xi,xj,r,dx,S;
  double S11,S12,S13,S22,S23,S33;
  char buf[256];  
  int i,nmax=1000;
  
  if (argc < 4) {
    fprintf(stderr,"Usage: %s xi xj r\n",argv[0]);
    exit(1);
  }
  xi = atof(argv[1]);
  xj = atof(argv[2]);
  dx = 0.001;
  r  = 0;
  for(i=0; (i<=nmax); i++) {
    S11 = Coulomb_SS(r,1,1,xi,xj);
    S12 = Coulomb_SS(r,1,2,xi,xj);
    S13 = Coulomb_SS(r,1,3,xi,xj);
    S22 = Coulomb_SS(r,2,2,xi,xj);
    S23 = Coulomb_SS(r,2,3,xi,xj);
    S33 = Coulomb_SS(r,3,3,xi,xj);
    printf("%8.3f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f\n",
	   r,S11,S12,S13,S22,S23,S33);
    r += dx;
  }
  return 0;
}

