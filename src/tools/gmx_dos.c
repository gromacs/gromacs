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
#include <stdio.h>
#include <math.h>

#include "confio.h"
#include "copyrite.h"
#include "gmx_fatal.h"
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
#include "correl.h"
#include "gmx_ana.h"
#include "gmx_fft.h"

int gmx_dos(int argc,char *argv[])
{
  const char *desc[] = {
    "[TT]g_dos[tt] computes the Density of States from a simulations.",
    "In order for this to be meaningful the velocities must be saved",
    "in the trajecotry with sufficiently high frequency such as to cover",
    "all vibrations. For flexible systems that would be around a few fs",
    "between saving. Properties based on the DoS are printed on the",
    "standard output."
  };
  const char *bugs[] = {
    "This program needs a lot of memory: total usage equals the number of atoms times 3 times number of frames times 4 (or 8 when run in double precision)."
  };
  FILE       *fp;
  t_topology top;
  int        ePBC=-1;
  t_trxframe fr;
  matrix     box;
  int        gnx;
  atom_id    *index;
  char       *grpname;
  char       title[256];
  real       t0,t1,m;
  t_trxstatus *status;
  int        teller,n_alloc,i,j,k,l,fftcode;
  double     dt;
  real       **c1,**dos,mi,beta;
  output_env_t oenv;
  gmx_fft_t  fft;
    
  static     gmx_bool bVerbose=TRUE;
  static     real Temp=298.15;
  t_pargs pa[] = {
    { "-v", FALSE, etBOOL, {&bVerbose},
      "Be Verbose" },
    { "-T", FALSE, etBOOL, {&Temp},
      "Temperature in the simulation" }
  };

  t_filenm  fnm[] = {
    { efTRN, "-f",    NULL,   ffREAD  },
    { efTPS, "-s",    NULL,   ffREAD }, 
    { efNDX, NULL,    NULL,   ffOPTRD },
    { efXVG, "-acf",  "acf",  ffWRITE },
    { efXVG, "-dos",  "dos",  ffWRITE }
  };
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,
		    asize(bugs),bugs,&oenv);

  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,NULL,NULL,box,
		TRUE);
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpname);
  gnx = gnx*DIM;
  
  /* Correlation stuff */
  snew(c1,gnx);
  for(i=0; (i<gnx); i++)
    c1[i]=NULL;
  
  read_first_frame(oenv,&status,ftp2fn(efTRN,NFILE,fnm),&fr,TRX_NEED_V);
  t0=fr.time;
      
  n_alloc=0;
  teller=0;
  do {
    if (teller >= n_alloc) {
      n_alloc+=100;
      for(i=0; i<gnx; i++)
	srenew(c1[i],n_alloc);
    }
    for(i=0; i<gnx; i+=DIM) {
      c1[i+XX][teller] = fr.v[index[i/DIM]][XX];
      c1[i+YY][teller] = fr.v[index[i/DIM]][YY];
      c1[i+ZZ][teller] = fr.v[index[i/DIM]][ZZ];
    }

    t1=fr.time;

    teller ++;
  } while (read_next_frame(oenv,status,&fr));
  
  close_trj(status);

  dt = (t1-t0)/(teller-1);
  if (bVerbose)
    printf("Going to do %d fourier transforms of length %d. Hang on.\n",
	   gnx,teller);
  do_autocorr(NULL,oenv,NULL,teller,gnx,c1,dt,eacNormal,FALSE);
  snew(dos,2);
  snew(dos[0],teller);
  snew(dos[1],teller);
  if (bVerbose)
    printf("Going to merge the ACFs into the mass-weighted total ACF\n");
  for(i=0; (i<gnx); i+=DIM) {
    mi = top.atoms.atom[i/DIM].m;
    for(j=0; (j<teller/2); j++) 
      dos[0][j] += mi * (c1[i+XX][j] + c1[i+YY][j] + c1[i+ZZ][j]);
  }
  fp = xvgropen(opt2fn("-acf",NFILE,fnm),"Mass-weighted velocity ACF",
		"Time (ps)","C(t)",oenv);
  for(j=0; (j<teller/2); j++) 
    fprintf(fp,"%10g  %10g\n",j*dt,dos[0][j]);
  fclose(fp);
    
  if ((fftcode = gmx_fft_init_1d_real(&fft,teller/2,GMX_FFT_FLAG_NONE)) != 0) 
    gmx_fatal(FARGS,"gmx_fft_init_1d_real returned %d",fftcode);
  if ((fftcode = gmx_fft_1d_real(fft,GMX_FFT_REAL_TO_COMPLEX,
				 (void *)dos[0],(void *)dos[1])) != 0)
    gmx_fatal(FARGS,"gmx_fft_1d_real returned %d",fftcode);

  beta = 1/(2*Temp*BOLTZ);
  fp = xvgropen(opt2fn("-dos",NFILE,fnm),"Density of states",
		"\\f{12}n\\f{4} (1/ps)","\\f{4}S(\\f{12}n\\f{4})",oenv);
  for(j=0; (j<teller/2); j++) 
    fprintf(fp,"%10g  %10g\n",(j == 0) ? 0 : (j*1.0/dt),dos[1][j*2]*beta);
  fclose(fp);
    
  do_view(oenv,ftp2fn(efXVG,NFILE,fnm),"-nxy");
  
  thanx(stderr);
  
  return 0;
}
