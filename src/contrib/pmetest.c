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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "copyrite.h"
#include "main.h"
#include "tpxio.h"
#include "statutil.h"
#include "futil.h"
#include "vec.h"
#include "ewald_util.h"

static void do_my_pme(FILE *fp,bool bVerbose,t_inputrec *ir,
		      rvec x[],rvec f[],real charge[],matrix box,
		      t_commrec *cr,t_nsborder *nsb,t_nrnb *nrnb,
		     real ewaldcoeff)
{
  real   ener;
  tensor vir;
  
  clear_mat(vir);
  
  ener = do_pme(fp,bVerbose,ir,x,f,charge,box,cr,
		nsb,nrnb,vir,ewaldcoeff,FALSE);
  fprintf(fp,"Energy: %12.5e\n",ener);
  pr_rvecs(fp,0,"virial",vir,DIM);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "The pmetest program tests the scaling of the PME code. When only given",
    "a tpr file it will compute PME for one frame. When given a trajectory",
    "it will do so for all the frames in the trajectory. Before the PME",
    "routine is called the coordinates are sorted along the X-axis."
  };
  t_commrec    *cr,*mcr;
  static t_filenm fnm[] = {
    { efTPX, NULL,      NULL,       ffREAD  },
    { efTRN, "-o",      NULL,       ffWRITE },
    { efLOG, "-g",      "pmetest",  ffWRITE },
    { efTRX, "-x",      NULL,       ffOPTRD }
  };
#define NFILE asize(fnm)

  /* Command line options ! */
  static bool bVerbose=FALSE;
  static bool bOptFFT=FALSE;
  static int  ewald_geometry=0;
  static int  nnodes=1;
  static int  nthreads=1;
  static int  pme_order=4;
  static rvec grid = { -1, -1, -1 };
  static real rc   = 0.9;
  static real dtol = 1e-4;
  static t_pargs pa[] = {
    { "-np",      FALSE, etINT, {&nnodes},
      "Number of nodes, must be the same as used for grompp" },
    { "-nt",      FALSE, etINT, {&nthreads},
      "Number of threads to start on each node" },
    { "-v",       FALSE, etBOOL,{&bVerbose},  
      "Be loud and noisy" },
    { "-grid",    FALSE, etRVEC,{&grid},
      "Number of grid cells in X, Y, Z dimension (if -1 use from tpr)" },
    { "-order",   FALSE, etINT, {&pme_order},
      "Order of the PME spreading algorithm" }
  };
  t_inputrec  ir;
  t_topology  top;
  t_tpxheader tpx;
  t_nrnb      nrnb;
  t_nsborder  *nsb;
  char        title[STRLEN];
  int         natoms,step,status,i;
  real        t,lambda,ewaldcoeff;
  rvec        *x,*f;
  real        *charge;
  matrix      box;
  
  cr = init_par(&argc,&argv);

  if (MASTER(cr))
    CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,
		    PCA_KEEP_ARGS | PCA_NOEXIT_ON_ARGS | PCA_BE_NICE |
		    PCA_CAN_SET_DEFFNM | (MASTER(cr) ? 0 : PCA_QUIET),
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  bVerbose = bVerbose && MASTER(cr);
    
#ifndef USE_MPI
  if (nnodes > 1) 
    gmx_fatal(FARGS,"GROMACS compiled without MPI support - can't do parallel runs");
#endif
#ifndef USE_THREADS
  if(nthreads > 1)
    gmx_fatal(FARGS,"GROMACS compiled without threads support - can only use one thread");
#endif

  open_log(ftp2fn(efLOG,NFILE,fnm),cr);

  /* Read tpr file etc. */
  read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&tpx,FALSE,NULL,NULL);
  snew(x,tpx.natoms);
  snew(f,tpx.natoms);
  snew(charge,tpx.natoms);
  read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,&ir,
	   box,&natoms,x,NULL,NULL,&top);
  /* Charges */
  snew(charge,tpx.natoms);
  for(i=0; (i<tpx.natoms); i++) 
    charge[i] = top.atoms.atom[i].q;
  
  /* Grid stuff */
  if (opt2parg_bSet("-grid",asize(pa),pa)) {
    ir.nkx = grid[XX];
    ir.nky = grid[YY];
    ir.nkz = grid[ZZ];
  }
  if ((ir.nkx <= 0) || (ir.nky <= 0) || (ir.nkz <= 0))
    gmx_fatal(FARGS,"PME grid = %d %d %d",ir.nkx,ir.nky,ir.nkz);
  
  init_pme(stdlog,cr,ir.nkx,ir.nky,ir.nkz,pme_order,
	   tpx.natoms,bOptFFT,ewald_geometry);
  init_nrnb(&nrnb);
  snew(nsb,1);
  /* Add parallellization code here */
  snew(nsb,1);
  
  ewaldcoeff = calc_ewaldcoeff(rc,dtol);
  fprintf(stdlog,"Ewald parameters:\n"
	  "rc         = %12.5e\n"
	  "dtol       = %12.5e\n"
	  "ewaldcoeff = %12.5e\n"
	  "grid       = %4d %4d %4d\n",
	  rc,dtol,ewaldcoeff,ir.nkx,ir.nky,ir.nkz);
	  
  fprintf(stdlog,"Results based on tpr file %s\n",ftp2fn(efTPX,NFILE,fnm));
  fflush(stdlog);
  
  do_my_pme(stdlog,bVerbose,&ir,x,f,charge,box,cr,nsb,&nrnb,ewaldcoeff);
  if (ftp2bSet(efTRX,NFILE,fnm)) {
    fprintf(stdlog,"Results based on trx file %s\n",ftp2fn(efTRX,NFILE,fnm));
    sfree(x);
    natoms = read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box); 
    if (natoms != top.atoms.nr)
      gmx_fatal(FARGS,"natoms in trx = %d, in tpr = %d",natoms,top.atoms.nr);
    do {
      do_my_pme(stdlog,bVerbose,&ir,x,f,charge,box,cr,nsb,&nrnb,ewaldcoeff);
    } while (read_next_x(status,&t,natoms,x,box));
    close_trx(status);
  }
  
  if (gmx_parallel_env)
    gmx_finalize(cr);

  if (MASTER(cr))
    thanx(stderr);
  
  return 0;
}

