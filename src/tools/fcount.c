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
 * GROningen MAchine for Chemical Simulation
 */
#include <stdio.h>
#include <math.h>
#include "string2.h"
#include "smalloc.h"
#include "statutil.h"

#define MAX 50000

static void my_usage(char *prog,char *arg)
{
  fprintf(stderr," Usage: %s [-p] [-s] [-f dist] [ -c #columns ]"
	  " naver < infile > outfile\n",prog);
  exit(1);
}

int main(int argc, char *argv[]) 
{
  double *x,**y,value=0.001;
  double *yav,yyy,yyy2,ymin,ymax,aver,var,sd,fdist;
  int    i,j,k,nav=100,ncol=1;
  char   buf[STRLEN];
  bool   bSilent=FALSE;
  bool   bFreq=FALSE;
  int    nn,*nfreq;
   
  fdist=0.0;
  for(i=1; (i<argc); i++) {
    if (argv[i][0] == '-')
      switch (argv[i][1]) {
      case 'p':
	value=1.0;
	break;
      case 'f':
	bFreq=TRUE;
	fdist=dscan(argc,argv,&i);
	break;
      case 'c':
	ncol=iscan(argc,argv,&i);
	break;
      case 's':
	bSilent=TRUE;
	break;
      default:
	usage(argv[0],argv[i]);
      }
    else
      nav=atoi(argv[i]);
  }
  if (!bSilent) {
    fprintf(stderr,"Will average stdin with %d columns, over %d points\n",
	    ncol,nav);
    if (bFreq)
      fprintf(stderr,"Will calc frequencies of occurance < %g\n",fdist);
  }
  snew(x,MAX);
  snew(y,ncol);
  snew(nfreq,ncol+1);
  for(i=0; (i<ncol); i++) {
    snew(y[i],MAX);
  }
  snew(yav,MAX);
  i=0;
  do {
    if (scanf("%s",buf) == 1) {
      if ((buf[0] != '@') && (buf[0] != '#')) {
	sscanf(buf,"%lf",&x[i]);
	nn=0;
	for(k=0; (k<ncol); k++) {
	  scanf("%lf",&(y[k][i]));
	  if (y[k][i] < fdist) {
	    nn++;
	  }
	}
	nfreq[nn]++;
	if (i >= nav) {
	  printf("%10e",x[i-nav/2]*value);
	  for(k=0; (k<ncol); k++) {
	    yav[k]=0;
	    for(j=i-nav+1; (j<=i); j++)
	      yav[k]+=y[k][j];
	    yav[k]/=nav;
	    printf("  %10e",yav[k]);
	  }
	  printf("\n");
	}
	i++;
      }
      else {
	while (getc(stdin) != '\n')
	  ;
      }
    }
    else
      break;
  } while ((int)strlen(buf) > 0);
  if (!bSilent)
    fprintf(stderr,"%3s  %12s  %12s  %12s  %12s  %12s\n",
	    "i","min","aver","max","var","sd");
  for(k=0; (k<ncol); k++) {
    aver=y[k][0];
    ymin=aver;
    ymax=aver;
    yyy2=aver*aver;
    for(j=1; (j<i); j++) {
      yyy=y[k][j];
      aver+=yyy;
      yyy2+=(yyy*yyy);
      if (yyy < ymin)
	ymin=yyy;
      else if (yyy > ymax)
	ymax=yyy;
    }
    aver/=i;
    var=yyy2/i-aver*aver;
    sd=sqrt(var);
    fprintf(stderr,"%3d  %12g  %12g  %12g  %12g  %12g\n",
	    k,ymin,aver,ymax,var,sd);
  }
  if (bFreq) {
    fprintf(stderr,"Number of Contacts < %g nm\n",fdist);
    for(k=0; (k<=ncol); k++)
      fprintf(stderr,"%2d: %5d times\n",k,nfreq[k]);
  }
  return 0;
}


