/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_g_analyze_c = "$Id$";

#include <math.h>
#include <string.h>
#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "fatal.h"
#include "vec.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "txtdump.h"
#include "gstat.h"
#include "xvgr.h"

real **read_val(char *fn,bool bHaveT,int nsets_in,
		int *nset,int *nval,real *t0,real *dt)
{
  FILE *fp;
  char line0[4096],*line,*format;
  int  a,narg,n,sin,set,nchar;
  double dbl,tend;
  bool bEndOfSet;
  real **val;

  val = NULL;
  fp  = ffopen(fn,"r");
  for(sin=0; sin<nsets_in; sin++) {
    if (nsets_in == 1)
      narg = 0;
    else 
      narg = bHaveT ? 2 : 1;
    n = 0;
    bEndOfSet = FALSE;
    while (!bEndOfSet && fgets(line0,STRLEN-1,fp)) {
      line = line0;
      bEndOfSet = (line[0] == '&');
      if ((line[0] != '#') && (line[0] != '@') && !bEndOfSet) {
	a = 0;
	while (((a<narg) || ((nsets_in==1) && (n==0))) && 
	       (line[0] != '\n') && sscanf(line,"%lf%n",&dbl,&nchar)) {
	  if (sin) {
	    if (!bHaveT || (a>0))
	      set = sin;
	    else
	      set = -1;
	  } else {
	    if (!bHaveT)
	      set = a;
	    else
	      set = a-1;
	  }
	  if (n==0) {
	    if (nsets_in == 1)
	      narg++;
	    if (set == -1)
	      *t0 = dbl;
	    else {
	      *nset = set+1;
	      srenew(val,*nset);
	      val[set] = NULL;
	    }
	  }
	  if (set == -1)
	    tend = dbl;
	  else {
	    if (n % 100 == 0)
	      srenew(val[set],n+100);
	    val[set][n] = (real)dbl;
	  }
	  a++;
	  line += nchar;
	}
	n++;
	if (a != narg)
	  fprintf(stderr,"Invalid line in %s: '%s'\n",fn,line0);
      }
    }
    if (sin==0) {
      *nval = n;
      if (!bHaveT)
	*dt = 1.0;
      else
	*dt = (real)(tend-*t0)/(n-1.0);
    } else {
      if (n < *nval) {
	fprintf(stderr,"Set %d is shorter (%d) than the previous set (%d)\n",
		sin+1,n,*nval);
	*nval = n;
	fprintf(stderr,"Will use only the first %d points of every set\n",
		*nval);
      }
    }
  }
  fclose(fp);
  
  return val;
}

void histogram(char *distfile,real binwidth,int n, int nset, real **val)
{
  FILE *fp;
  int  i,s;
  real min,max;
  int  nbin;
  real *histo;

  min=val[0][0];
  max=val[0][0];
  for(s=0; s<nset; s++)
    for(i=0; i<n; i++)
      if (val[s][i] < min)
	min = val[s][i];
      else if (val[s][i] > max)
	max = val[s][i];
  
  if (-min > max)
    max = -min;
  nbin = (int)(max/binwidth)+1;
  fprintf(stderr,"Making distributions with %d bins\n",2*nbin+1);
  snew(histo,2*nbin+1);
  fp = xvgropen(distfile,"Distribution","","");
  for(s=0; s<nset; s++) {
    for(i=0; i<2*nbin+1; i++)
      histo[i] = 0;
    for(i=0; i<n; i++)
      histo[nbin+(int)(floor(val[s][i]/binwidth+0.5))]++;
    for(i=0; i<2*nbin+1; i++)
      fprintf(fp," %g  %g\n",(i-nbin)*binwidth,(real)histo[i]/(n*binwidth));
    if (s<nset-1)
      fprintf(fp,"&\n");
  }
  fclose(fp);
}

void average(char *avfile,char **avbar_opt,
	     int n, int nset,real **val,real t0,real dt)
{
  FILE *fp;
  int  i,s;
  real av,var,err;
  char c;
  
  c = avbar_opt[0][0];

  fp = ffopen(avfile,"w");
  if ((c == 'e') && (nset == 1))
    c = 'n';
  if (c != 'n') 
    fprintf(fp,"@TYPE xydy\n");
  
  
  for(i=0; i<n; i++) {
    av = 0;
    for(s=0; s<nset; s++)
      av += val[s][i];
    av /= nset;
    fprintf(fp," %g %g",t0+dt*i,av);
    var = 0;
    if (c != 'n') {
      for(s=0; s<nset; s++)
	var += sqr(val[s][i]-av);
      if (c == 's')
	err = sqrt(var/nset);
      else
	err = sqrt(var/(nset*(nset-1)));
      fprintf(fp," %g",err);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

void estimate_error(char *eefile,int resol,int n,int nset,
		    real *av,real **val,real dt)
{
  FILE *fp;
  int log2max,rlog2,bs,prev_bs,nb;
  int s,i,j;
  real blav,var;
  char **leg;

  log2max = (int)(log(n)/log(2));

  snew(leg,nset);
  for(s=0; s<nset; s++) {
    snew(leg[s],STRLEN);
    sprintf(leg[s],"av %f",av[s]);
  }

  fp = xvgropen(eefile,"Error estimates","Block size (time)","Error estimate");
  fprintf(fp,
	  "@ subtitle \"using block averaging, total time %g (%d points)\"\n",
	  n*dt,n);
  xvgr_legend(fp,nset,leg);
  for(s=0; s<nset; s++)
    sfree(leg[s]);
  sfree(leg);

  for(s=0; s<nset; s++) {
    prev_bs = 0;
    for(rlog2=resol*log2max; rlog2>=2*resol; rlog2--) {
      bs = n*pow(0.5,(real)rlog2/(real)resol);
      if (bs != prev_bs) {
	nb = 0;
	i = 0;
	var = 0;
	while (i+bs <= n) {
	  blav=0;
	  for (j=0; j<bs; j++) {
	    blav += val[s][i];
	  i++;
	  }
	  var += sqr(av[s] - blav/bs);
	  nb++;
	}
	fprintf(fp," %g %g\n",bs*dt,sqrt(var/(nb*(nb-1))));
      }
      prev_bs = bs;
    }
    if (s < nset)
      fprintf(fp,"&\n");
  }

  fclose(fp);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_analyze reads an ascii file and analyzes data sets.",
    "A line in the input file may start with a time",
    "(see option [TT]-time[tt]) and any number of y values may follow.",
    "Multiple sets can also be",
    "read when they are seperated by & (option [TT]-n[tt]),",
    "in this case only one y value is read from each line.",
    "All lines starting with # and @ are skipped.",
    "All analyses can also be done for the derivative of a set",
    "(option [TT]-d[tt]).[PAR]",
    "Option [TT]-ac[tt] produces the autocorrelation function(s).[PAR]",
    "Option [TT]-msd[tt] produces the mean square displacement(s).[PAR]",
    "Option [TT]-dist[tt] produces distribution plot(s).[PAR]",
    "Option [TT]-av[tt] produces the average over the sets,",
    "optionally with error bars ([TT]-errbar[tt]).[PAR]",
    "Option [TT]-ee[tt] produces error estimates using block averaging.",
    "A set is divided in a number of blocks and averages are calculated for",
    "each block. The error for the total average is calculated from the",
    "variance between the block averages. These errors are plotted as a",
    "function of the block size. For a good error estimate the block size",
    "should be at least as large as the correlation time, but possibly much",
    "larger.[PAR]"
  };
  static real frac=0.5,binwidth=0.1;
  static bool bHaveT=TRUE,bDer=FALSE,bSubAv=FALSE,bAverCorr=FALSE;
  static int  nsets_in=1,d=1,resol=8;

  static char *avbar_opt[] = { NULL, "none", "stddev", "error", NULL };

  t_pargs pa[] = {
    { "-time", FALSE, etBOOL, {&bHaveT},
      "Expect a time in the input" },
    { "-n", FALSE, etINT, {&nsets_in},
      "Read # sets seperated by &" },
    { "-d", FALSE, etBOOL, {&bDer},
	"Use the derivative" },
    { "-dp",  FALSE, etINT, {&d}, 
      "HIDDENThe derivative is the difference over # points" },
    { "-bw", FALSE, etREAL, {&binwidth},
      "Binwidth for the distribution" },
    { "-errbar", FALSE, etENUM, {&avbar_opt},
      "Error bars for the average" },
    { "-resol", FALSE, etINT, {&resol},
      "HIDDENResolution for the block averaging, block size increases with"
    " a factor 2^(1/#)" },
    { "-subav", FALSE, etBOOL, {&bSubAv},
      "Subtract the average before autocorrelating" },
    { "-oneacf", FALSE, etBOOL, {&bAverCorr},
      "Calculate one ACF over all sets" }
  };
#define NPA asize(pa)

  FILE     *out;
  int      n,nlast,s,nset,i,t;
  real     **val,t0,dt,tot,*av;
  char     *acfile,*msdfile,*distfile,*avfile,*eefile;
  
  t_filenm fnm[] = { 
    { efXVG, "-f",    "graph",    ffREAD   },
    { efXVG, "-ac",   "autocorr", ffOPTWR  },
    { efXVG, "-msd",  "msd",      ffOPTWR  },
    { efXVG, "-dist", "distr",    ffOPTWR  },
    { efXVG, "-av",   "average",  ffOPTWR  },
    { efXVG, "-ee",   "errest",   ffOPTWR  }
  }; 
#define NFILE asize(fnm) 

  int     npargs;
  t_pargs *ppa;

  npargs = asize(pa); 
  ppa    = add_acf_pargs(&npargs,pa);
  
  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,0,TRUE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL); 

  acfile   = opt2fn_null("-ac",NFILE,fnm);
  msdfile  = opt2fn_null("-msd",NFILE,fnm);
  distfile = opt2fn_null("-dist",NFILE,fnm);
  avfile   = opt2fn_null("-av",NFILE,fnm);
  eefile   = opt2fn_null("-ee",NFILE,fnm);
  if (!acfile && !msdfile && !distfile && !avfile && !eefile) {
    fprintf(stderr,"Please use one of the output file options\n");
    exit(0);
  }

  val=read_val(opt2fn("-f",NFILE,fnm),bHaveT,nsets_in,&nset,&n,&t0,&dt);
  fprintf(stderr,"Read %d sets of %d points, dt = %g\n",nset,n,dt);
  if (bDer) {
    fprintf(stderr,"Calculating the derivative as (f[i+%d]-f[i])/(%d*dt)\n",
	    d,d);
    n -= d;
    for(s=0; s<nset; s++)
      for(i=0; i<n; i++)
	val[s][i] = (val[s][i+d]-val[s][i])/(d*dt);
  }	

  snew(av,nset);
  for(s=0; s<nset; s++) {
    for(i=0; i<n; i++)
      av[s] += val[s][i];
    av[s] /= n;
    fprintf(stderr,"Average of set %d: %g\n",s+1,av[s]); 
  }

  if (msdfile) {
    out=xvgropen(msdfile,"Mean square displacement",
		 "time (ps)","MSD (nm\\S2\\N)");
    nlast = (int)(n*frac);
    for(s=0; s<nset; s++) {
      for(t=0; t<=nlast; t++) {
	if (t % 100 == 0)
	  fprintf(stderr,"\r%d",t);
	tot=0;
	for(i=0; i<n-t; i++)
	  tot += sqr(val[s][i]-val[s][i+t]); 
	tot /= (real)(n-t);
	fprintf(out," %g %8g\n",dt*t,tot);
      }
      if (s<nset-1)
	fprintf(out,"&\n");
    }
    fclose(out);
    fprintf(stderr,"\r%d, time=%g\n",t-1,(t-1)*dt);
  }
  
  if (distfile)
    histogram(distfile,binwidth,n,nset,val);
  
  if (avfile)
    average(avfile,avbar_opt,n,nset,val,t0,dt);

  if (eefile)
    estimate_error(eefile,resol,n,nset,av,val,dt);

  if (acfile) {
    if (bSubAv) 
      for(s=0; s<nset; s++)
	for(i=0; i<n; i++)
	  val[s][i] -= av[s];
    do_autocorr(acfile,"Autocorrelation",n,nset,val,dt,
		eacNormal,bAverCorr,NULL,NULL);
  }

  return 0;
}
  
