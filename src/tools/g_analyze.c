#include <stdio.h>
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
    while (fgets(line0,STRLEN-1,fp) && !bEndOfSet) {
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
    "optionally with error bars ([TT]-errbar[tt]).",
  };
  static real frac=0.5,binwidth=0.1;
  static bool bHaveT=TRUE,bDer=FALSE,bAver=FALSE;
  static int  nsets_in=1,d=1;

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
    { "-aver", FALSE, etBOOL, {&bAver},
      "Average all autocorrelation functions" },
    { "-bw", FALSE, etREAL, {&binwidth},
      "binwidth for the distribution" },
    { "-errbar", FALSE, etENUM, {&avbar_opt},
      "Error bars for the average" }
  };
#define NPA asize(pa)

  FILE     *out;
  int      n,nlast,s,nset,i,t;
  real     **val,t0,dt,tot;
  char     *acfile,*msdfile,*distfile,*avfile;
  
  t_filenm fnm[] = { 
    { efXVG, "-f",    "graph",    ffREAD   },
    { efXVG, "-ac",   "autocorr", ffWRITE  },
    { efXVG, "-msd",  "msd",      ffOPTWR  },
    { efXVG, "-dist", "distr",    ffOPTWR  },
    { efXVG, "-av",   "average",  ffOPTWR  }
  }; 
#define NFILE asize(fnm) 

  int     npargs;
  t_pargs *ppa;

  npargs = asize(pa); 
  ppa    = add_acf_pargs(&npargs,pa);
  
  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,0,TRUE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL); 

  acfile = opt2fn_null("-ac",NFILE,fnm);
  msdfile = opt2fn_null("-msd",NFILE,fnm);
  distfile = opt2fn_null("-dist",NFILE,fnm);
  avfile = opt2fn_null("-av",NFILE,fnm);
  if (!acfile && !msdfile && !distfile && !avfile) {
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
  
  if (acfile)
    do_autocorr(acfile,"Autocorrelation",n,nset,val,dt,
		eacNormal,bAver,NULL,NULL);

  return 0;
}
  
