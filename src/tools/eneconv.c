/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
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
static char *SRCID_eneconv_c = "$Id$";

#include <string.h>
#include <math.h>
#include "string2.h"
#include "typedefs.h"
#include "smalloc.h"
#include "statutil.h"
#include "disre.h"
#include "names.h"
#include "copyrite.h"
#include "macros.h"
#include "enxio.h"

typedef struct {
  real      t,tstart;
  t_energy  *e;
  t_drblock *dr;
} t_fixe;

bool bRmod(double a,double b)
{
  int iq;
  double tol = 1e-6;
  
  iq = ((1.0+tol)*a)/b;
  
  if (fabs(a-b*iq) <= tol*a)
    return TRUE;
  else
    return FALSE;
}

static void set_fixe(t_fixe *fixe, real t, real tstart,
		     int nre, t_energy ee[], t_drblock *dr)
{
  fixe->t      = t;
  fixe->tstart = tstart;
  snew(fixe->e,nre);
  memcpy(fixe->e,ee,sizeof(ee[0])*nre);
  fixe->dr=dr;
}

int cmpfix(const void *a,const void *b)
{
  real ts,tt;

  /* First sort on starting time of the simulation */  
  ts = ((t_fixe *) a)->tstart - ((t_fixe *) b)->tstart;
  tt = ((t_fixe *) a)->t      - ((t_fixe *) b)->t;
  
  if (ts == 0) {
    if (tt < 0.0)
      return -1;
    else if (tt > 0.0)
      return 1;
    else
      return 0;
  }
  else {
    if (ts < 0.0)
      return -1;
    else if (ts > 0.0)
      return 1;
    else
      return 0;
  }
}

bool bRgt(double a,double b)
{
  double tol = 1e-6;
  
  if ( a > (b - tol*(a+b)) )
    return TRUE;
  else
    return FALSE;
}

void analyse_energies(int nre,int nrf,t_fixe f[])
{
  real     tstart,t;
  t_energy *e0,*elast;
  int      i,j,nrkill;
  
  if (!nrf)
    return;
    
  tstart = f[0].tstart;
  
  /* First throw out redundant (double) energies */
  nrkill=0;
  for(i=0; (i<nrf); i++) {
    t = f[i].t;
    if (t != -1) {
      /* Check whether this is data from the same simulation  */
      if (f[i].tstart != tstart) {
	tstart = f[i].tstart;
	/* Walk back until we have the same time */
	while ((i > 0) && bRgt(f[i-1].t,t) ) {
	  /* Set time to -1 to indicate that this should not be stored */
	  f[i-1].t = -1;
	  nrkill ++;
	  i--;
	}
	/* Now either i == 0, or f[i-1].t <= t */
      }
    }
  }
  fprintf(stderr,"There are %d double energy entries out of %d\n",
	  nrkill,nrf);
  
  /* Now set the averages etc. */
  snew(e0,nre);
  snew(elast,nre);
  tstart = -1;
  for(i=0; (i<nrf); i++) {
    t = f[i].t;
    if (t != -1) {
      memcpy(elast,f[i].e,nre*sizeof(elast[0]));
      if ((tstart == -1) || (f[i].tstart != tstart)) {
	tstart = f[i].tstart;
	memcpy(e0,f[i].e,nre*sizeof(e0[0]));
      }
      
    }    
  }
  sfree(e0);
  sfree(elast);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "with only [TT]-o[tt] specified:[BR]",
    "repairs energy files, by checking the contents and",
    "concatenating the files in the right order. In case of double", 
    "time frames one of these is skipped.",
    "The input files are taken from the command line, such that the command",
    "[TT]fix_ene -o fixed.ene *.ene[tt] should do the trick.",
    "(note that in this case the -o option should come first and ",
    "should include a filename)",
    "[PAR]",
    "With [TT]-f[tt] and [TT]-od[tt] specified: converts edr to ene.",
    "[PAR]",
    "With [TT]-fd[tt] and [TT]-o[tt] specified: converts ene to edr."
  };
  
  static char *bugs[] = {
    "Will use rediculous amounts of memory for large files, which "
    "slows things down."
  };

  int       in,out;
  t_energy  *ee=NULL;
  int       step,nrt,nre,nresav;
  real      t,tsav,tstart;
  bool      bSetTime;
  char      **fnms,fn[STRLEN];
  int       nfile;
  t_fixe    *fixe=NULL;
  int       *set,i,j;
  char      **enm=NULL;
  t_drblock *dr=NULL;

  t_filenm fnm[] = {
    { efENX, "-f", NULL,    ffREAD  },
    { efENX, "-o", "fixed", ffOPTWR },
  };
#define NFILE asize(fnm)
  
  static real t0=-1, timestep=-1, delta_t=0.0, toffset=0;
  t_pargs pa[] = {
    { "-dt",       FALSE,  etREAL, &delta_t,
      "only write out frame when (t MOD delta_t) == offset" },
    { "-offset",   FALSE, etREAL, &toffset,
      "time offset for -dt option" },
    { "-t0",       FALSE, etREAL, &t0, 
      "change starting time" },
    { "-timestep", FALSE, etREAL, &timestep, 
      "change timestep between frames" }
  };
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_NOEXIT_ON_ARGS,FALSE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);
  
  snew(fnms,argc-1);
  nfile=0;
  bSetTime=FALSE;
  for(i=1; (i<argc); i++)
    if (strcmp(fn,argv[i]) != 0)
      fnms[nfile++]=argv[i];
  if (nfile) {
    fprintf(stderr,"Files to fix:");
    for (i=0; (i<nfile); i++)
      fprintf(stderr,"%s ",fnms[i]);
    fprintf(stderr,"\n");
  }
  bSetTime=((t0 != -1) || (timestep != -1));
  
  if (!opt2bSet("-f",NFILE,fnm)) {
  
    strcpy(fn,opt2fn("-o",NFILE,fnm));
    out=open_enx(fn,"w");
    
    nresav=0;
    nrt=0;
    
    if (nfile==0)
      fatal_error(0,"No input files to fix!");
    
    for(i=0; (i<nfile); i++) {
      in = open_enx(fnms[i],"r");
      fprintf(stderr,"Opened %s\n",fnms[i]);
      do_enxnms(in,&nre,&enm);
      tstart=-1;
      if (i == 0) {
	do_enxnms(out,&nre,&enm);
	nresav=nre;
	snew(ee,nre);
      }
      else if (nre != nresav)
	fatal_error(0,"Energy files don't match, different number"
		    " of energies (%s)",fnms[i]);
      step=0;
      while (do_enx(in,&t,&step,&nre,ee,dr)) {
	fprintf(stderr,"\rRead frame %d, time %g",nrt,t);
	if (tstart == -1)
	  tstart = t;
	srenew(fixe,++nrt);
	set_fixe(&fixe[nrt-1],t,tstart,nre,ee,dr);
      }
      close_enx(in);
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"Sorting energies...\n");
    qsort(fixe,nrt,sizeof(fixe[0]),cmpfix);
    fprintf(stderr,"Analysing energies...\n");
    analyse_energies(nre,nrt,fixe);
    
    tsav=-1;
    j=0;
    if ((t0==-1) && (timestep>0))
      t0=tstart;
    fprintf(stderr,"Now writing %s...\n",fn);
    for(i=0; (i<nrt); i++) { 
      if (bSetTime) {
	if (timestep>0)
	  fixe[i].t=t0+j*timestep;
	else if (t0>=0) 
	  fixe[i].t += (t0-tstart);
      }
      if ( ( bRmod(fixe[i].t-toffset,delta_t) || (delta_t==0) ) && 
	   ( bSetTime || ( fixe[i].t != tsav ) ) ) {
	do_enx(out,&(fixe[i].t),&i,&nre,fixe[i].e,fixe[i].dr);
	fprintf(stderr,"\rWrite frame %d, t %f",j,fixe[i].t);
	j++;
      }
      tsav=fixe[i].t;
    }
    close_enx(out);
    fprintf(stderr,"\n");
  } else {
    in = open_enx(opt2fn("-f",NFILE,fnm),"r");
    
    dr=get_drblock();
    do_enxnms(in,&nre,&enm);
    snew(ee,nre);
    
    out = open_enx(opt2fn("-o",NFILE,fnm),"w");
    do_enxnms(out,&nre,&enm);
    
    tstart=-1;
    i=0;
    while (do_enx(in,&t,&j,&nre,ee,dr))  {
      if (tstart==-1) {
	tstart=t;
	if ((t0==-1) && (timestep>0))
	  t0=tstart;
      }
      if ((j%10) == 0)
	fprintf(stderr,"\rRead frame %d, time %g",j,t);
      if (timestep>0)
	t=t0+j*timestep;
      else if (t0>=0) 
	t += (t0-tstart);
      if ( bRmod(t-toffset,delta_t) || (delta_t==0) ) {
	fprintf(stderr,"  ->  write %d, time %g",i,t);
	do_enx(out,&t,&j,&nre,ee,dr);
	i++;
      }
      j++;
    }
    fprintf(stderr,"\n");
    
  };
  
  thanx(stdout);
  
  return 0;
}
