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
 * Grunge ROck MAChoS
 */
static char *SRCID_g_velacc_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "bondf.h"
#include "confio.h"
#include "copyrite.h"
#include "fatal.h"
#include "futil.h"
#include "gstat.h"
#include "macros.h"
#include "maths.h"
#include "physics.h"
#include "rdgroup.h"
#include "smalloc.h"
#include "statutil.h"
#include "string.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "typedefs.h"
#include "vec.h"
#include "strdb.h"
#include "xvgr.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_velacc computes the velocity autocorrelation function"
  };
  static int  nf = 10;
  t_pargs pa[] = {
    { "-nframes",   FALSE, etINT,  &nf,
      "Number of frames in your trajectory" }
  };

  FILE       *fp;
  rvec       *v;
  int        gnx;
  atom_id    *index;
  char       *grpname;
  char       title[256];
  real       t,t0,t1,lambda;
  matrix     box;
  int        status,natoms,teller,i,tel3;
  real       **c1;
  
#define NHISTO 360
    
  t_filenm  fnm[] = {
    { efTRJ, "-f",    NULL,   ffREAD },
    { efXVG, "-o",    "vac",  ffWRITE },
    { efNDX, NULL,    NULL,   ffREAD }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL);

  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);
  sprintf(title,"Velocity Autocorrelation Function for %s",grpname);
  
  /* Correlation stuff */
  snew(c1,gnx);
  if (nf < 0) 
    fatal_error(0,"No frames (%d) in trajectory ? DIY!\n",nf);
  fprintf(stderr,"Going to malloc %d bytes!\n",gnx*DIM*nf*sizeof(c1[0][0]));
  for(i=0; (i<gnx); i++)
    snew(c1[i],DIM*nf);
  
  natoms=read_first_v(&status,opt2fn("-f",NFILE,fnm),&t,&v,box);
  t0=t;
      
  teller=0;
  do {
    if ((teller % 10) == 0)
      fprintf(stderr,"\rt=%.2f",t);
    
    if (teller >= nf)
      break;
    tel3=3*teller;
    for(i=0; (i<gnx); i++) {
      c1[i][tel3+XX]=v[index[i]][XX];
      c1[i][tel3+YY]=v[index[i]][YY];
      c1[i][tel3+ZZ]=v[index[i]][ZZ];
    }
    
    teller ++;
  } while (read_next_v(status,&t,natoms,v,box));
  fprintf(stderr,"\n");
  close_trj(status);

  t1=t;
  
  do_autocorr(ftp2fn(efXVG,NFILE,fnm),"Velocity Autocorrelation Function",
	      teller,gnx,c1,(t1-t0)/teller,eacVector,
	      TRUE,NULL,NULL);
  
  xvgr_file(ftp2fn(efXVG,NFILE,fnm),"-nxy");
  
  thanx(stdout);
  
  return 0;
}
