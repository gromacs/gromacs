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
#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "enxio.h"
#include "statutil.h"
#include "names.h"
#include "copyrite.h"
#include "macros.h"

typedef struct {
  real     t,tstart;
  t_energy *e;
} t_fixe;

static void set_fixe(t_fixe *fixe,real t,real tstart,
		     int nre,t_energy ee[])
{
  fixe->t      = t;
  fixe->tstart = tstart;
  snew(fixe->e,nre);
  memcpy(fixe->e,ee,sizeof(ee[0])*nre);
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
	while ((i > 0) && (f[i-1].t > t)) {
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
    "fix_ene repairs energy files, by checking the contents and",
    "concatenating the files in the right order. In case of double", 
    "time frames one of these is skipped.[BR]",
    "The input files are taken from the command line, such that the command",
    "[TT]fix_ene -o fixed.ene *.ene[tt] should do the trick."
  };

  int      in,out;
  t_energy *ee;
  int      step,nrt,nre,nresav;
  real     t,tsav,t0,dt,tstart;
  bool     bSetTime;
  char     **fnms;
  int      nfile;
  t_fixe   *fixe=NULL;
  int      *set,i,j;
  char     **enm;
  char     *outfile;
  t_filenm fnm[] = {
    { efENX, "-o", "fixed", ffWRITE }
  };
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_NOEXIT_ON_ARGS,
		    FALSE,1,fnm,0,NULL,asize(desc),desc,0,NULL);

  outfile = ftp2fn(efENE,1,fnm);
  out = open_enx(outfile,"w");
  
  nresav=0;
  nrt=0;
  
  snew(fnms,argc-1);
  nfile=0;
  bSetTime=FALSE;
  for(i=1; (i<argc); i++) {
    if (strcmp(argv[i],"-t") == 0) {
      t0=dscan(argc,argv,&i);
      dt=dscan(argc,argv,&i);
      bSetTime=TRUE;
    }
    else if (strcmp(outfile,argv[i]) == 0) {
      fprintf(stderr,"Skipping %s\n",outfile);
      continue;
    }
    else {
      fnms[nfile++]=argv[i];
    }
  }
  if (bSetTime && (nfile != 1)) {
    fatal_error(0,"Can only set the time for one file at a time!");
  }
  
  for(i=0; (i<nfile); i++) {
    in = open_enx(fnms[i],"r");
    fprintf(stderr,"Opened %s\n",fnms[i]);
    do_enxnms(in,&nre,&enm);
    tstart=-1;
    if (nresav == 0) {
      do_enxnms(out,&nre,&enm);
      nresav=nre;
      snew(ee,nre);
    }
    else if (nre != nresav)
      fatal_error(0,"Energy files don't match, different number"
		  " of energies");
    while (do_enx(in,&t,&step,&nre,ee,NULL)) {
      fprintf(stderr,"\rt: %.2f",t);
      if (tstart == -1)
	tstart = t;
      srenew(fixe,++nrt);
      set_fixe(&fixe[nrt-1],t,tstart,nre,ee);
    }
    close_enx(in);
    fprintf(stderr,"\n");
  }
  fprintf(stderr,"Analysing energies...\n");
  qsort(fixe,nrt,sizeof(fixe[0]),cmpfix);
  analyse_energies(nre,nrt,fixe);
  
  tsav=-1;
  j=0;
  fprintf(stderr,"Now writing fixed.ene...\n");
  for(i=0; (i<nrt); i++) { 
    if (bSetTime) {
      fixe[i].t=t0+j*dt;
      do_enx(out,&fixe[i].t,&i,&nre,fixe[i].e,NULL);
      fprintf(stderr,"\rframe: %d",j);
      j++;
    }
    else if (fixe[i].t != tsav) {
      do_enx(out,&fixe[i].t,&i,&nre,fixe[i].e,NULL);
      fprintf(stderr,"\rframe: %d",j);
      j++;
    }
    tsav=fixe[i].t;
  }
  close_enx(out);
  fprintf(stderr,"\n");
  
  thanx(stdout);
  
  return 0;
}
