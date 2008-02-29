/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include "string2.h"
#include "smalloc.h"
#include "statutil.h"
#include "vec.h"

static void my_usage(char *prog,char *arg)
{
  fprintf(stderr," Usage: %s [-p] [-s] [-c #columns]"
	  " naver < infile > outfile\n",prog);
  fprintf(stderr,"-p means picoseconds rather than nanoseconds\n");
  fprintf(stderr,"-s means silent rather than verbose\n");
  fprintf(stderr,"Don't ever use argument %s again!\n",arg);
  exit(1);
}

void lsq_y_ax_b_double(int n, double x[], double y[], double *a, double *b)
{
  int    i;
  double yx,xx,sx,sy;

  yx=xx=sx=sy=0.0;
  for (i=0; (i < n); i++) {
    yx+=y[i]*x[i];
    xx+=x[i]*x[i];
    sx+=x[i];
    sy+=y[i];
  }
  *a=(n*yx-sy*sx)/(n*xx-sx*sx);
  *b=(sy-(*a)*sx)/n;
}

int main(int argc, char *argv[]) 
{
  double *x,**y,value=0.001;
  double *yav,yyy,yyy2,ymin,ymax,aver,var,sd;
  int    i,j,k,nav=100,ncol=1,MAX=50000;
  char   buf[STRLEN];
  bool   bSilent=FALSE,bVerySilent=FALSE;
  double   lsq_a,lsq_b,rms_res;
  
  for(i=1; (i<argc); i++) {
    if (argv[i][0] == '-')
      switch (argv[i][1]) {
      case 'p':
	value=1.0;
	break;
      case 'm':
	MAX=iscan(argc,argv,&i);
	break;
      case 'c':
	ncol=iscan(argc,argv,&i);
	break;
      case 's':
	bSilent=TRUE;
	break;
      case 'S':
	bVerySilent=bSilent=TRUE;
	break;
      default:
	my_usage(argv[0],argv[i]);
      }
    else
      nav=atoi(argv[i]);
  }
  if (!bSilent)
    fprintf(stderr,"Will average stdin with %d columns, over %d points\n",
	    ncol,nav);
  
  
  snew(x,MAX);
  snew(y,ncol);
  for(i=0; (i<ncol); i++) {
    snew(y[i],MAX);
  }
  snew(yav,MAX);
  i=0;
  do {
    if (scanf("%s",buf) == 1) {
      if ((buf[0] != '@') && (buf[0] != '#')) {
	sscanf(buf,"%lf",&x[i]);
	for(k=0; (k<ncol); k++)
	  scanf("%lf",&(y[k][i]));
	if (i >= nav) {
	  if (!bVerySilent) 
	    printf("%10e",x[i-nav/2]*value);
	  for(k=0; (k<ncol); k++) {
	    yav[k]=0;
	    for(j=i-nav+1; (j<=i); j++)
	      yav[k]+=y[k][j];
	    yav[k]/=nav;
	    if (!bVerySilent) 
	      printf("  %10e",yav[k]);
	  }
	  if (!bVerySilent) 
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
  } while (((int)strlen(buf) > 0) && (i<MAX));


  if (!bSilent)
    fprintf(stderr,"%3s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n",
	    "i","min","aver","max","var","sd","drift","fluc");
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

    lsq_y_ax_b_double(i,x,y[k],&lsq_a,&lsq_b);
    rms_res=0;
    for(j=0; (j<i);j++)
      rms_res+=sqr(y[k][j]-(lsq_a*x[j]+lsq_b));
    rms_res=sqrt(rms_res/i);

    fprintf(stderr,"%3d  %12g  %12g  %12g  %12g  %12g  %12g  %12g\n",
	    k,ymin,aver,ymax,var,sd,lsq_a,rms_res);
  }
  return 0;
}


