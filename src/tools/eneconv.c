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
  bool      bSave;
  int       sim_nr;
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
		     int nre, t_energy ee[], t_drblock *dr,int nr)
{
  fixe->t      = t;
  fixe->tstart = tstart;
  fixe->bSave  = TRUE;
  fixe->sim_nr = nr;
  snew(fixe->e,nre);
  memcpy(fixe->e,ee,sizeof(ee[0])*nre);
  fixe->dr=dr;
}

static bool same_time(real t1,real t2)
{
  const real tol=1e-5;

  return (fabs(t1-t2) < tol);
}

int cmpfix(const void *a,const void *b)
{
  real ta,tb;

  ta = ((t_fixe *) a)->t;
  tb = ((t_fixe *) b)->t;
  
  if (same_time(ta,tb))
    return 0;
  else if (ta < tb)
    return -1;
  else 
    return 1;
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
  real     t;
  t_energy *e0,*elast;
  int      i,i0,j,nrkill;
  
  if (!nrf)
    return;
    
  /* First throw out redundant (double) energies. 
   * They are sorted on time, so that is pretty easy.
   */
  nrkill=0;
  for(i=0; (i<nrf-1); ) {
    i0 = i;
    while (((i+1) < nrf) && same_time(f[i0].t,f[i+1].t)) {
      f[i+1].bSave = FALSE;
      nrkill ++;
      i++;
    }
    if (i0 == i)
      i++;
  }
  fprintf(stderr,"There are %d double energy entries out of %d\n",
	  nrkill,nrf);
  fprintf(stderr,"Not setting averages (not implemented)\n");
  /* Now set the averages etc. 
  snew(e0,nre);
  snew(elast,nre);
  sim_start = -1;
  for(i=0; (i<nrf); i++) {
    t = f[i].t;
    if (t != -1) {
      memcpy(elast,f[i].e,nre*sizeof(elast[0]));
      if ((sim_start == -1) || (f[i].sim_nr != sim_start)) {
	sim_start = f[i].sim_nr;
	memcpy(e0,f[i].e,nre*sizeof(e0[0]));
      }
      
    }    
  }
  sfree(e0);
  sfree(elast);*/
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "With [TT]-f[tt] is [IT]not[it] specified:[BR]",
    "repairs energy files, by checking the contents and concatenating",
    "the files in the right order. In case of double time frames one of",
    "these is skipped. The input files are taken from the command line,",
    "such that the command [TT]eneconv -o fixed.ene *.ene[tt] should do",
    "the trick. (note that in this case the -o option must come first and ",
    "must include a filename) [PAR]",
    "With [TT]-f[tt] specified:[BR]",
    "reads one energy file and writes another, applying the [TT]-dt[tt],",
    "[TT]-offset[tt], [TT]-t0[tt] and [TT]-timestep[tt] options and",
    "converting to a different format if necessary (indicated by file",
    "extentions)."
  };
  
  static char *bugs[] = {
    "Will use ridiculous amounts of memory for large files, which "
    "slows things down."
  };

  int       in,out;
  t_energy  *ee=NULL;
  int       step,nrt,nre,nresav;
  real      t,tsav,tstart;
  bool      bSetTime,bWriteThis;
  char      **fnms,fn[STRLEN];
  int       nfile;
  t_fixe    *fixe=NULL;
  int       *set,i,j,ndr;
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
      "only write out frame when t MOD delta_t = offset" },
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
  /* If dt or starting time is specified, we have to change it */
  bSetTime=((t0 != -1) || (timestep != -1));
  
  if (!opt2bSet("-f",NFILE,fnm)) {
  
    if (nfile==0)
      fatal_error(0,"No input files to fix!");
      
    strcpy(fn,opt2fn("-o",NFILE,fnm));
    out=open_enx(fn,"w");
    
    nresav=0;
    nrt=0;
    
    for(i=0; (i<nfile); i++) {
      /* Loop over input energy files */
      in = open_enx(fnms[i],"r");
      fprintf(stderr,"Opened %s\n",fnms[i]);
      do_enxnms(in,&nre,&enm);
      tstart=-1;
      if (i == 0) {
	do_enxnms(out,&nre,&enm);
	nresav = nre;
	snew(ee,nre);
      }
      else if (nre != nresav)
	fatal_error(0,"Energy files don't match, different number"
		    " of energies (%s)",fnms[i]);
		    
      /* Now really start reading */
      step = 0;
      while (do_enx(in,&t,&step,&nre,ee,&ndr,dr)) {
	fprintf(stderr,"\rRead frame %d, time %g",nrt,t);
	/* Set the starting time for this energy file */
	if (tstart == -1)
	  tstart = t;
	srenew(fixe,++nrt);
	set_fixe(&fixe[nrt-1],t,tstart,nre,ee,dr,i);
      }
      close_enx(in);
      fprintf(stderr,"\n");
    }
    /* Now the fixe array contains all energies from all files */
    fprintf(stderr,"Sorting energies...\n");
    qsort(fixe,nrt,sizeof(fixe[0]),cmpfix);
    fprintf(stderr,"Analysing energies...\n");
    analyse_energies(nre,nrt,fixe);
    
    j=0;
    if ((t0==-1) && (timestep>0))
      t0=tstart;
    fprintf(stderr,"Now writing %s...\n",fn);
    for(i=0; (i<nrt); i++) { 
      /* Check whether we have to write this step */
      bWriteThis = (delta_t==0) || (bRmod(fixe[i].t-toffset,delta_t));

      if (bWriteThis && fixe[i].bSave) {
      
	if (bSetTime) {
	  if (timestep>0)
	    fixe[i].t  = t0+j*timestep;
	  else if (t0>=0) 
	    fixe[i].t += (t0-tstart);
	}
	do_enx(out,&(fixe[i].t),&i,&nre,fixe[i].e,&ndr,fixe[i].dr);
	fprintf(stderr,"\rWrite frame %d, t %f",j,fixe[i].t);
	j++;
      }
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
    while (do_enx(in,&t,&j,&nre,ee,&ndr,dr))  {
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
      if ((delta_t == 0) || (bRmod(t-toffset,delta_t))) {
	if ((j%10) == 0)
	  fprintf(stderr,"  ->  write %d, time %g",i,t);
	do_enx(out,&t,&j,&nre,ee,&ndr,dr);
	i++;
      }
      j++;
    }
    fprintf(stderr,"\n");
    
  };
  
  thanx(stdout);
  
  return 0;
}
