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
 * GROwing Monsters And Cloning Shrimps
 */
#include <math.h>
#include <stdio.h>
#include <time.h>

#include "typedefs.h"
#include "macros.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/xvgr.h"
#include "copyrite.h"
#include "mdrun.h"
#include "main.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fft/fft.h"
#include "gromacs/math/gmxcomplex.h"

#ifdef GMX_MPI
#include "gromacs/fft/parallel_3dfft.h"
#endif

#include "fftgrid.h"


int main(int argc,char *argv[])
{
  int       mmm[] = { 8, 10, 12, 15, 16, 18, 20, 24, 25, 27, 30, 32, 36, 40,
		      45, 48, 50, 54, 60, 64, 72, 75, 80, 81, 90, 100 };
  int       nnn[] = { 24, 32, 48, 60, 72, 84, 96 };
#define NNN asize(nnn)
  FILE      *fp,*fplog;
  int       *niter;
  int       i,j,n,nit,ntot,n3,rsize;
  double    t,nflop,start;
  double    *rt,*ct;
  t_fftgrid *g;
  t_commrec *cr;
  static gmx_bool bReproducible = FALSE;
  static int  nnode    = 1;
  static int  nitfac  = 1;
  t_pargs pa[] = {
    { "-reproducible",   FALSE, etBOOL, {&bReproducible}, 
      "Request binary reproducible results" },
    { "-np",    FALSE, etINT, {&nnode},
      "Number of NODEs" },
    { "-itfac", FALSE, etINT, {&nitfac},
      "Multiply number of iterations by this" }
  };
  static t_filenm fnm[] = {
    { efLOG, "-g", "fft",      ffWRITE },
    { efXVG, "-o", "fft",      ffWRITE }
  };
#define NFILE asize(fnm)
  
  cr = init_par(&argc,&argv);
  if (MASTER(cr))
    CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv, PCA_CAN_SET_DEFFNM,
		    NFILE,fnm,asize(pa),pa,0,NULL,0,NULL);
  gmx_log_open(ftp2fn(efLOG,NFILE,fnm),cr,1,0,&fplog);

  snew(niter,NNN);
  snew(ct,NNN);
  snew(rt,NNN);
  rsize = sizeof(real);
  for(i=0; (i<NNN); i++) {
    n  = nnn[i];
    if (n < 16)
      niter[i] = 50;
    else if (n < 26)
      niter[i] = 20;
    else if (n < 51)
      niter[i] = 10;
    else
      niter[i] = 5;
    niter[i] *= nitfac;
    nit = niter[i];
    
    if (MASTER(cr))
      fprintf(stderr,"\r3D FFT (%s precision) %3d^3, niter %3d     ",
	      (rsize == 8) ? "Double" : "Single",n,nit);
    
    g  = mk_fftgrid(n,n,n,NULL,NULL,cr,bReproducible);

    if (PAR(cr))
      start = time(NULL);
    else
      start_time();
    for(j=0; (j<nit); j++) {
      gmxfft3D(g,GMX_FFT_REAL_TO_COMPLEX,cr);
      gmxfft3D(g,GMX_FFT_COMPLEX_TO_REAL,cr);
    }
    if (PAR(cr)) 
      rt[i] = time(NULL)-start;
    else {
      update_time();
      rt[i] = node_time();
    }
    done_fftgrid(g);
    sfree(g);
  }
  if (MASTER(cr)) {
    fprintf(stderr,"\n");
    fp=xvgropen(ftp2fn(efXVG,NFILE,fnm),
		"FFT timings","n^3","t (s)");
    for(i=0; (i<NNN); i++) {
      n3 = 2*niter[i]*nnn[i]*nnn[i]*nnn[i];
      fprintf(fp,"%10d  %10g\n",nnn[i],rt[i]/(2*niter[i]));
    }
    gmx_fio_fclose(fp);
  }
  return 0;
}
