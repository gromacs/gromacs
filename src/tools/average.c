/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
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
 * Gnomes, ROck Monsters And Chili Sauce
 */
static char *SRCID_average_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "string2.h"
#include "smalloc.h"
#include "statutil.h"
#include "vec.h"

#define MAX 50000

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
  int    i,j,k,nav=100,ncol=1;
  char   buf[STRLEN];
  bool   bSilent=FALSE;
  double   lsq_a,lsq_b,rms_res;
  
  for(i=1; (i<argc); i++) {
    if (argv[i][0] == '-')
      switch (argv[i][1]) {
      case 'p':
	value=1.0;
	break;
      case 'c':
	ncol=iscan(argc,argv,&i);
	break;
      case 's':
	bSilent=TRUE;
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


