/*
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

#include <math.h>
#include <string.h>
#include "sysstuff.h"
#include "physics.h"
#include "typedefs.h"
#include "smalloc.h"
#include "futil.h"
#include "statutil.h"
#include "copyrite.h"
#include "index.h"
#include "macros.h"
#include "gmx_fatal.h"
#include "xvgr.h"
#include "gstat.h"
#include "vec.h"
#include "gmx_ana.h"


int gmx_rotacf(int argc,char *argv[])
{
  const char *desc[] = {
    "g_rotacf calculates the rotational correlation function",
    "for molecules. Three atoms (i,j,k) must be given in the index",
    "file, defining two vectors ij and jk. The rotational acf",
    "is calculated as the autocorrelation function of the vector",
    "n = ij x jk, i.e. the cross product of the two vectors.",
    "Since three atoms span a plane, the order of the three atoms",
    "does not matter. Optionally, controlled by the -d switch, you can",
    "calculate the rotational correlation function for linear molecules",
    "by specifying two atoms (i,j) in the index file.",
    "[PAR]",
    "EXAMPLES[PAR]",
    "g_rotacf -P 1 -nparm 2 -fft -n index -o rotacf-x-P1",
    "-fa expfit-x-P1 -beginfit 2.5 -endfit 20.0[PAR]",
    "This will calculate the rotational correlation function using a first",
    "order Legendre polynomial of the angle of a vector defined by the index",
    "file. The correlation function will be fitted from 2.5 ps till 20.0 ps",
    "to a two parameter exponential.",


    ""
  };
  static gmx_bool bVec    = FALSE,bAver=TRUE;

  t_pargs pa[] = {
    { "-d",   FALSE, etBOOL, {&bVec},
      "Use index doublets (vectors) for correlation function instead of triplets (planes)" },
    { "-aver",FALSE, etBOOL, {&bAver},
      "Average over molecules" }
  };

  t_trxstatus *status;
  int        isize;
  atom_id    *index;
  char       *grpname;
  rvec       *x,*x_s;
  matrix     box;
  real       **c1;
  rvec       xij,xjk,n;
  int        i,m,teller,n_alloc,natoms,nvec,ai,aj,ak;
  unsigned long mode;
  real       t,t0,t1,dt;
  gmx_rmpbc_t  gpbc=NULL;
  t_topology *top;
  int        ePBC;
  t_filenm   fnm[] = {
    { efTRX, "-f", NULL,  ffREAD  },
    { efTPX, NULL, NULL,  ffREAD },
    { efNDX, NULL, NULL,  ffREAD  },
    { efXVG, "-o", "rotacf",  ffWRITE }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  output_env_t oenv;
  
  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL,&oenv);
  
  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&isize,&index,&grpname);
  
  if (bVec) 
    nvec = isize/2;
  else
    nvec = isize/3;
  
  if (((isize % 3) != 0) && !bVec)
    gmx_fatal(FARGS,"number of index elements not multiple of 3, "
		"these can not be atom triplets\n");
  if (((isize % 2) != 0) && bVec)
    gmx_fatal(FARGS,"number of index elements not multiple of 2, "
		"these can not be atom doublets\n");
  
  top=read_top(ftp2fn(efTPX,NFILE,fnm),&ePBC);
  
  snew(c1,nvec);
  for (i=0; (i<nvec); i++)
    c1[i]=NULL;
  n_alloc=0;

  natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  snew(x_s,natoms);
  
  gpbc = gmx_rmpbc_init(&(top->idef),ePBC,natoms,box);

  /* Start the loop over frames */
  t1 = t0 = t;
  teller  = 0;
  do {
    if (teller >= n_alloc) {
      n_alloc+=100;
      for (i=0; (i<nvec); i++)
	srenew(c1[i],DIM*n_alloc);
    }
    t1 = t;
    
    /* Remove periodicity */
    gmx_rmpbc_copy(gpbc,natoms,box,x,x_s);
  
    /* Compute crossproducts for all vectors, if triplets.
     * else, just get the vectors in case of doublets.
     */
    if (bVec == FALSE) {
      for (i=0; (i<nvec); i++) {
	ai=index[3*i];
	aj=index[3*i+1];
	ak=index[3*i+2];
	rvec_sub(x_s[ai],x_s[aj],xij);
	rvec_sub(x_s[aj],x_s[ak],xjk);
	cprod(xij,xjk,n);
	for(m=0; (m<DIM); m++)
	  c1[i][DIM*teller+m]=n[m];
      }
    }
    else {
      for (i=0; (i<nvec); i++) {
	ai=index[2*i];
	aj=index[2*i+1];
	rvec_sub(x_s[ai],x_s[aj],n);
	for(m=0; (m<DIM); m++)
	  c1[i][DIM*teller+m]=n[m];
      }
    }
    /* Increment loop counter */
    teller++;
  } while (read_next_x(oenv,status,&t,natoms,x,box));  
  close_trj(status); 
  fprintf(stderr,"\nDone with trajectory\n");
  
  gmx_rmpbc_done(gpbc);


  /* Autocorrelation function */
  if (teller < 2)
    fprintf(stderr,"Not enough frames for correlation function\n");
  else {
    dt=(t1 - t0)/(teller-1);
    
    mode = eacVector;
    
    do_autocorr(ftp2fn(efXVG,NFILE,fnm),oenv,"Rotational Correlation Function",
		teller,nvec,c1,dt,mode,bAver);
  }

  do_view(oenv,ftp2fn(efXVG,NFILE,fnm),NULL);
    
  thanx(stderr);
    
  return 0;
}
