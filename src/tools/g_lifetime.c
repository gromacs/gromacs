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
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_g_lifetime_c = "$Id$";

#include "typedefs.h"
#include "string2.h"
#include "copyrite.h"
#include "smalloc.h"
#include "futil.h"
#include "xvgr.h"
#include "macros.h"
#include "statutil.h"

#define MAXLINE 10000

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_lifetime reads a hydrogen-bond frequency matrix as computed",
    "by g_hbond, and computes the autocorrelation function.",
    "This should give an estimate for the lifetime of hydrogen bonds."
  };
  static int  nobonds=-1;
  static real dt=1;
  t_pargs pa [] = {
    { "-n", FALSE, etINT, &nobonds,
      "number of hydrogen bonds considered, if not given (-1) all hydrogen bonds are taken into account." },
    { "-dt",FALSE, etREAL,&dt,
      "time between frames (ps)" }
  };
  static char *bugs[] = {
    "This is really S L O W"
  };

  FILE *fp,*fps; /* file pointers */
  char line[MAXLINE];    /* read frame */
  int  *n;       /* binary 1 0 array */
  real  *ac;      /* autocorrelation function */
  int i,j;       /* just counters */
  int t,tau;     /* time and tau counter */
  int nbonds;    /* the number of hydrogen bonds (cols ) */
  int nframes;   /* the number of time frames ( rows ) */
  int bond,nlast;
  real t2;
  t_filenm fnm[] = {
    efDAT, "-f", "freq", FALSE,
    efXVG, "-o", "lifetime", FALSE
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);
  
  /* open file */
  fp = ftp2FILE(efDAT,NFILE,fnm,"r");

  /* read number of rows and columns */
  fgets2(line,MAXLINE,fp);
  sscanf(line,"%d %d",&nbonds,&nframes);
  if (nbonds > MAXLINE) 
    fatal_error(0,"Too many hbonds (%d), change MAXLINE in %s.c",
		nbonds,argv[0]);
		
  fprintf(stderr,"nframes = %5d\nnbonds = %5d\n",nframes,nbonds);
  if ((nobonds == -1) || (nobonds > nbonds)) 
    nobonds=nbonds;
  
  fprintf(stderr,"I will use %d hbonds out of %d\n",nobonds,nbonds);
  
  /* malloc read frame and data arrays */
  snew(n,nframes+1);
  snew(ac,nframes+1);

  /* read and process the data */
  for(bond=0;(bond<nobonds);bond++) {
    fprintf(stderr,"\rbond no %d of %d           ",bond+1,nbonds);
    /* the current hydrogen bond is i */
    rewind(fp);
    fgets2(line,MAXLINE,fp);
    /* read all frames of i */
    for(t=0;(t<nframes);t++) {
      fgets2(line,MAXLINE,fp);
      n[t]=line[bond]-'0';
    }

    /* determine contribution */
    for(t=0;(t<nframes);t++) {
      for(tau=0;(tau<nframes-t);tau++) {
	ac[tau]+=(real)(n[t]*n[t+tau]);
      }
    }

  }

  /* close file */
  fclose(fp);

  /* print autocorrelation */
  fp=xvgropen(ftp2fn(efXVG,NFILE,fnm),"Lifetime of HBonds",
	      "Time (ps)","C(t)");
  for(t=1;(t<nframes);t++)
    if (ac[t] < 0.5*ac[0]) {
      t2=(t-1+(0.5*ac[0]-ac[t-1])/(ac[t]-ac[t-1]))*dt;
      fprintf(fp,"@ subtitle \"t\\s1/2\\N = %.3f ps\"\n",t2);
      break;
    }

  for(nlast=nframes; (t>0); t--)
    if (ac[nlast-1] != 0)
      break;
      
  for(t=0;(t<nlast);t++)
    fprintf(fp,"%10e   %10.5e\n",t*dt,ac[t]/ac[0]);
  fclose(fp);
  
  xvgr_file(ftp2fn(efXVG,NFILE,fnm),NULL);
  
  thanx(stdout);
  
  return 0;
}





