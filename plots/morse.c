#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void main(int argc,char *argv[])
{
  int i;
  double x,y,Dij,kij,beta,bij;
  
  if (argc == 4) {
    Dij = atof(argv[1]);
    kij = atof(argv[2]);
    bij = atof(argv[3]);
  }
  else {
    fprintf(stderr,"Usage %s Dij kij bij\n",argv[0]);
    exit(1);
  }
  beta = sqrt(kij/(2*Dij));
  for(i=0; (i<=100); i++) {
    x = (i*1.0)/100;
    y = Dij*pow((1-exp(-beta*(x-bij))),2.0);
    printf("%10g  %10g\n",x,y);
  }
}
