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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_g_velacc_c = "$Id$";

#include <stdio.h>
#include <math.h>
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

  t_trxframe fr;
  int        gnx;
  atom_id    *index;
  char       *grpname;
  char       title[256];
  real       t0,t1;
  int        status,teller,n_alloc,i,tel3;
  real       **c1;
  
#define NHISTO 360
    
  t_filenm  fnm[] = {
    { efTRN, "-f",    NULL,   ffREAD },
    { efNDX, NULL,    NULL,   ffREAD },
    { efXVG, "-o",    "vac",  ffWRITE }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  CopyRight(stderr,argv[0]);
  npargs = 0;
  ppa    = add_acf_pargs(&npargs,NULL);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL);

  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);
  sprintf(title,"Velocity Autocorrelation Function for %s",grpname);
  
  /* Correlation stuff */
  snew(c1,gnx);
  for(i=0; (i<gnx); i++)
    c1[i]=NULL;
  
  read_first_frame(&status,ftp2fn(efTRN,NFILE,fnm),&fr,TRX_NEED_V);
  t0=fr.time;
      
  n_alloc=0;
  teller=0;
  do {
    if (teller >= n_alloc) {
      n_alloc+=100;
      for(i=0; (i<gnx); i++)
	srenew(c1[i],DIM*n_alloc);
    }
    tel3=3*teller;
    for(i=0; (i<gnx); i++) {
      c1[i][tel3+XX]=fr.v[index[i]][XX];
      c1[i][tel3+YY]=fr.v[index[i]][YY];
      c1[i][tel3+ZZ]=fr.v[index[i]][ZZ];
    }
    
    teller ++;
  } while (read_next_frame(status,&fr));
  close_trj(status);

  t1=fr.time;

  do_autocorr(ftp2fn(efXVG,NFILE,fnm),"Velocity Autocorrelation Function",
	      teller,gnx,c1,(t1-t0)/(teller-1),eacVector,TRUE);
  
  do_view(ftp2fn(efXVG,NFILE,fnm),"-nxy");
  
  thanx(stderr);
  
  return 0;
}
