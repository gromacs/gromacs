/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#include <stdio.h>

#define MAX 10000

int main(int argc, char *argv[]) 
{
  double x[MAX],y_1[MAX],y_2[MAX],y_3[MAX],value=0.001;
  double yav_1,yav_2,yav_3;
  int    i,j,nav;
  char   buf[257];
 
  if (argc<2) {
    fprintf(stderr," Usage: %s [-ps] naver < infile > outfile\n",argv[0]);
    exit(1);
  } 
   
  if ((argv[1][0]=='-')&&(argv[1][1]=='p')&&(argv[1][2]='s')) {
    value=1.0;
    nav=atoi(argv[2]);
  }
  else 
    nav=atoi(argv[1]);

  
  i=0;
  while (fgets(buf,256,stdin) != NULL) {
    if (buf[0] != '@') {
      sscanf(buf,"%lf %lf %lf %lf",&x[i],&y_1[i],&y_2[i],&y_3[i]);
      if (i >= nav) {
	yav_1=0;
	yav_2=0;
	yav_3=0;
	for(j=i-nav+1; (j<=i); j++) {
	  yav_1+=y_1[j];
	  yav_2+=y_2[j];
	  yav_3+=y_3[j];
	}
	yav_1/=nav;
	yav_2/=nav;
	yav_3/=nav;
	printf("%12g  %12g %12g %12g\n",x[i]*value,yav_1,yav_2,yav_3);
      }
      i++;
    }
  }
}


