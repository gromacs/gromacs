/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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
#include "index.h"
#include "smalloc.h"
#include "statutil.h"
#include "string.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "typedefs.h"
#include "vec.h"
#include "strdb.h"
#include "xvgr.h"

int gmx_velacc(int argc,char *argv[])
{
  static char *desc[] = {
    "g_velacc computes the velocity autocorrelation function.",
    "When the [TT]-s[tt] option is used, the momentum autocorrelation",
    "function is calculated.[PAR]",
    "With option [TT]-mol[tt] the momentum autocorrelation function of",
    "molecules is calculated. In this case the index group should consist",
    "of molecule numbers instead of atom numbers."
  };
  
  static bool bMol=FALSE;
  t_pargs pa[] = {
    { "-mol", FALSE, etBOOL, {&bMol},
      "Calculate vac of molecules" }
  };

  t_topology top;
  t_trxframe fr;
  matrix     box;
  bool       bTPS,bTop=FALSE;
  int        gnx;
  atom_id    *index,*a=NULL,*atndx=NULL,at;
  char       *grpname;
  char       title[256];
  real       t0,t1,m;
  int        status,teller,n_alloc,i,j,tel3;
  rvec       mv_mol;
  real       **c1;
  
#define NHISTO 360
    
  t_filenm  fnm[] = {
    { efTRN, "-f",    NULL,   ffREAD  },
    { efTPS, NULL,    NULL,   ffOPTRD }, 
    { efNDX, NULL,    NULL,   ffOPTRD },
    { efXVG, "-o",    "vac",  ffWRITE }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL);

  bTPS = bMol || ftp2bSet(efTPS,NFILE,fnm) || !ftp2bSet(efNDX,NFILE,fnm);

  if (bTPS) {
    bTop=read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,NULL,NULL,box,TRUE);
    get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);
  } else
    rd_index(ftp2fn(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);

  if (bMol) {
    if (!bTop)
      gmx_fatal(FARGS,"Need a topology to determine the molecules");
    a     = top.blocks[ebMOLS].a;
    atndx = top.blocks[ebMOLS].index;
  }
  
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
      for(i=0; i<gnx; i++)
	srenew(c1[i],DIM*n_alloc);
    }
    tel3=3*teller;
    if (bMol)
      for(i=0; i<gnx; i++) {
	clear_rvec(mv_mol);
	for(j=0; j<atndx[index[i]+1] - atndx[index[i]]; j++) {
	  at = a[atndx[index[i]]+j];
	  m  = top.atoms.atom[at].m;
	  mv_mol[XX] += m*fr.v[at][XX];
	  mv_mol[YY] += m*fr.v[at][YY];
	  mv_mol[ZZ] += m*fr.v[at][ZZ];
	}
	c1[i][tel3+XX]=mv_mol[XX];
	c1[i][tel3+YY]=mv_mol[YY];
	c1[i][tel3+ZZ]=mv_mol[ZZ];
      }
    else
      for(i=0; i<gnx; i++) {
	if (bTPS)
	  m = top.atoms.atom[index[i]].m;
	else
	  m = 1;
	c1[i][tel3+XX]=m*fr.v[index[i]][XX];
	c1[i][tel3+YY]=m*fr.v[index[i]][YY];
	c1[i][tel3+ZZ]=m*fr.v[index[i]][ZZ];
      }

    t1=fr.time;

    teller ++;
  } while (read_next_frame(status,&fr));
  close_trj(status);

  do_autocorr(ftp2fn(efXVG,NFILE,fnm),"Velocity Autocorrelation Function",
	      teller,gnx,c1,(t1-t0)/(teller-1),eacVector,TRUE);
  
  do_view(ftp2fn(efXVG,NFILE,fnm),"-nxy");
  
  thanx(stderr);
  
  return 0;
}
