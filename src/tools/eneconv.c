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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_eneconv_c = "$Id$";

#include <string.h>
#include "string2.h"
#include "typedefs.h"
#include "smalloc.h"
#include "enerio.h"
#include "statutil.h"
#include "disre.h"
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
    "with only [TT]-o[tt] specified: ",
    "repairs energy files, by checking the contents and",
    "concatenating the files in the right order. In case of double", 
    "time frames one of these is skipped.",
    "The input files are taken from the command line, such that the command",
    "[TT]fix_ene -o fixed.ene *.ene[tt] should do the trick.",
    "[PAR]",
    "With [TT]-f[tt] and [TT]-od[tt] specified: converts edr to ene.",
    "[PAR]",
    "With [TT]-fd[tt] and [TT]-o[tt] specified: converts ene to edr."
  };

  FILE     *in,*out;
  t_energy *ee;
  int      step,nrt,nre,nresav;
  real     t,tsav,tstart;
  bool     bSetTime,bF,bFD,bO,bOD,bFix,bEne2Edr,bEdr2Ene;
  char     **fnms,fn[STRLEN];
  int      nfile;
  t_fixe   *fixe=NULL;
  int      *set,i,j;
  char     **enm;
  t_filenm fnm[] = {
    { efENE, "-f", NULL,    ffOPTRD },
    { efEDR, "-fd", NULL,   ffOPTRD },
    { efENE, "-o", "fixed", ffOPTWR },
    { efEDR, "-od","fixed", ffOPTWR }
  };
#define NFILE asize(fnm)
  static real t0=-1, dt=-1;
  t_pargs pa[] = {
    { "-t",  FALSE, etREAL, &t0, "starting time" },
    { "-dt", FALSE, etREAL, &dt, "timestep" }
  };
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_NOEXIT_ON_ARGS,
		    FALSE,NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  
  bSetTime=((t0 != -1) || (dt != -1));
  bF =opt2bSet("-f",NFILE,fnm);
  bFD=opt2bSet("-fd",NFILE,fnm);
  bO =opt2bSet("-o",NFILE,fnm);
  bOD=opt2bSet("-od",NFILE,fnm);
  if (bF && bFD) fatal_error(0,"conflicting arguments: -f -fd");
  if (bO && bOD) fatal_error(0,"conflicting arguments: -o -od");
  bFix=bO && !( bOD || bF || bFD );
  bEne2Edr=!bFix && bF && bOD;
  bEdr2Ene=( !bFix && !bEne2Edr && bFD && bO );
  if (bFix) {
  
    strcpy(fn,opt2fn("-o",NFILE,fnm));
    out=ffopen(fn,"w");
    
    nresav=0;
    nrt=0;
    
    snew(fnms,argc-1);
    nfile=0;
    bSetTime=FALSE;
    for(i=1; (i<argc); i++) {
      if (strcmp(fn,argv[i]) == 0) {
	fprintf(stderr,"Skipping %s\n",fn);
	continue;
      }
      else {
	fnms[nfile++]=argv[i];
      }
    }
    if (nfile==0)
      fatal_error(0,"No input files to fix!");
    if (bSetTime && (nfile != 1))
      fatal_error(0,"Can only set the time for one file at a time!");
    
    for(i=0; (i<nfile); i++) {
      in=ffopen(fnms[i],"r");
      fprintf(stderr,"Opened %s\n",fnms[i]);
      rd_ener_nms(in,&nre,&enm);
      tstart=-1;
      if (nresav == 0) {
	wr_ener_nms(out,nre,enm);
	nresav=nre;
	snew(ee,nre);
      }
      else if (nre != nresav)
	fatal_error(0,"Energy files don't match, different number"
		    " of energies");
      while (rd_ener(in,&t,&step,ee,NULL)) {
	fprintf(stderr,"\rRead frame %d, time %g",nrt,t);
	if (tstart == -1)
	  tstart = t;
	srenew(fixe,++nrt);
	set_fixe(&fixe[nrt-1],t,tstart,nre,ee);
      }
      fclose(in);
      fprintf(stderr,"\n");
    }
    fprintf(stderr,"Analysing energies...\n");
    qsort(fixe,nrt,sizeof(fixe[0]),cmpfix);
    analyse_energies(nre,nrt,fixe);
    
    tsav=-1;
    j=0;
    if ((t0==-1) && (dt>0))
      t0=tstart;
    fprintf(stderr,"Now writing %s...\n",fn);
    for(i=0; (i<nrt); i++) { 
      if (bSetTime) {
	if (dt>0)
	  fixe[i].t=t0+j*dt;
	else if (t0>=0) 
	  fixe[i].t += (t0-tstart);
      }
      if ( bSetTime || (fixe[i].t != tsav) ) {
	wr_ener(out,fixe[i].t,i,nre,fixe[i].e,NULL);
	fprintf(stderr,"\rWrite frame %d, t %f",j,fixe[i].t);
	j++;
      }
      tsav=fixe[i].t;
    }
    fclose(out);
    fprintf(stderr,"\n");
    
  } else if (bEne2Edr) {
    FILE      *in;
    t_energy  *ee;
    t_drblock *dr;
    XDR       xdr;
    int       nre,j;
    real      t;
    char      **enm;
    
    in=ftp2FILE(efENE,NFILE,fnm,"r");
    
    dr=get_drblock();
    rd_ener_nms(in,&nre,&enm);
    snew(ee,nre);
    xdropen(&xdr,ftp2fn(efEDR,NFILE,fnm),"w");
    edr_nms(&xdr,&nre,&enm);
    
    tstart=-1;
    while (rd_ener(in,&t,&j,ee,dr))  {
      if (tstart==-1) {
	tstart=t;
	if ((t0==-1) && (dt>0))
	  t0=tstart;
      }
      if ((j%10) == 0)
	fprintf(stderr,"\rRead frame %d, time %g",j,t);
      if (dt>0)
	t=t0+j*dt;
      else if (t0>=0) 
	t += (t0-tstart);
      edr_io(&xdr,&t,&j,&nre,ee,dr);
      j++;
    }
    fprintf(stderr,"\n");
    
  } else if (bEdr2Ene) {
    FILE      *out;
    t_energy  *ee;
    t_drblock *dr;
    XDR       xdr;
    int       nre,j;
    real      t;
    char      **enm;
    
    tstart=-1;
    xdropen(&xdr,ftp2fn(efEDR,NFILE,fnm),"r");
    enm = NULL;
    edr_nms(&xdr,&nre,&enm);
    
    dr=get_drblock();
    snew(ee,nre);
    
    out=ftp2FILE(efENE,NFILE,fnm,"w");
    wr_ener_nms(out,nre,enm);
    
    j=0;
    while ( edr_io(&xdr,&t,&j,&nre,ee,dr)) {
      if (tstart==-1) {
	tstart=t;
	if ((t0==-1) && (dt>0))
	  t0=tstart;
      }
      if ((j%10) == 0)
	fprintf(stderr,"\rRead frame %d, time %g",j,t);
      if (dt>0)
	t=t0+j*dt;
      else if (t0>=0) 
	t += (t0-tstart);
      wr_ener(out,t,j,nre,ee,dr);
      j++;
    }
    fprintf(stderr,"\n");
    
  } else 
    fatal_error(0,"Nothing to be done....");
  
  thanx(stdout);
  
  return 0;
}
