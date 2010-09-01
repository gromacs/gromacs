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

#include "statutil.h"
#include "smalloc.h"
#include "macros.h"
#include "confio.h"
#include "copyrite.h"
#include "xvgr.h"
#include "index.h"
#include "do_fit.h"
#include "gmx_ana.h"
#include "gmx_fatal.h"


static real dointerp(int n,rvec x1[],rvec x2[],rvec xx[],
		    int I,int N,real first,real last)
{
  int    i,j;
  double fac,fac0,fac1;
  
  fac  = first + (I*(last-first))/(N-1);
  fac0 = 1-fac;
  fac1 = fac;
  for(i=0; (i<n); i++) 
    for(j=0; (j<DIM); j++)
      xx[i][j] = fac0*x1[i][j] + fac1*x2[i][j];
      
  return fac;
}

int gmx_morph(int argc,char *argv[])
{
  const char *desc[] = {
    "g_morph does a linear interpolation of conformations in order to",
    "create intermediates. Of course these are completely unphysical, but",
    "that you may try to justify yourself. Output is in the form of a ",
    "generic trajectory. The number of intermediates can be controlled with",
    "the -ninterm flag. The first and last flag correspond to the way of",
    "interpolating: 0 corresponds to input structure 1 while",
    "1 corresponds to input structure 2.",
    "If you specify first < 0 or last > 1 extrapolation will be",
    "on the path from input structure x1 to x2. In general the coordinates",
    "of the intermediate x(i) out of N total intermidates correspond to:[PAR]",
    "x(i) = x1 + (first+(i/(N-1))*(last-first))*(x2-x1)[PAR]",
    "Finally the RMSD with respect to both input structures can be computed",
    "if explicitly selected (-or option). In that case an index file may be",
    "read to select what group RMS is computed from."
  };
  t_filenm fnm[] = {
    { efSTX, "-f1", "conf1",  ffREAD },
    { efSTX, "-f2", "conf2",  ffREAD },
    { efTRX, "-o",  "interm", ffWRITE },
    { efXVG, "-or", "rms-interm", ffOPTWR },
    { efNDX, "-n",  "index",  ffOPTRD }
  };
#define NFILE asize(fnm)
  static  int  ninterm = 11;
  static  real first   = 0.0;
  static  real last    = 1.0;
  static  gmx_bool bFit    = TRUE;
  t_pargs pa [] = {
    { "-ninterm", FALSE, etINT,  {&ninterm},
      "Number of intermediates" },
    { "-first",   FALSE, etREAL, {&first},
      "Corresponds to first generated structure (0 is input x0, see above)" },
    { "-last",    FALSE, etREAL, {&last},
      "Corresponds to last generated structure (1 is input x1, see above)" },
    { "-fit",     FALSE, etBOOL, {&bFit},
      "Do a least squares fit of the second to the first structure before interpolating" }
  };
  const char *leg[] = { "Ref = 1\\Sst\\N conf", "Ref = 2\\Snd\\N conf" };
  FILE     *fp=NULL;
  int      i,isize,is_lsq,nat1,nat2;
  t_trxstatus *status;
  atom_id  *index,*index_lsq,*index_all,*dummy;
  t_atoms  atoms;
  rvec     *x1,*x2,*xx,*v;
  matrix   box;
  real     rms1,rms2,fac,*mass;
  char     title[STRLEN],*grpname;
  gmx_bool     bRMS;
  output_env_t oenv;
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    0,NULL,&oenv);
  get_stx_coordnum (opt2fn("-f1",NFILE,fnm),&nat1);
  get_stx_coordnum (opt2fn("-f2",NFILE,fnm),&nat2);
  if (nat1 != nat2)
    gmx_fatal(FARGS,"Number of atoms in first structure is %d, in second %d",
		nat1,nat2);
  
  init_t_atoms(&atoms,nat1,TRUE);
  snew(x1,nat1);
  snew(x2,nat1);
  snew(xx,nat1);
  snew(v,nat1);
  
  read_stx_conf(opt2fn("-f1",NFILE,fnm),title,&atoms,x1,v,NULL,box);
  read_stx_conf(opt2fn("-f2",NFILE,fnm),title,&atoms,x2,v,NULL,box);

  snew(mass,nat1);
  snew(index_all,nat1);
  for(i=0; (i<nat1); i++) {
    mass[i] = 1;
    index_all[i] = i;
  }
  if (bFit) {
    printf("Select group for LSQ superposition:\n");
    get_index(&atoms,opt2fn_null("-n",NFILE,fnm),1,&is_lsq,&index_lsq,
	      &grpname);
    reset_x(is_lsq,index_lsq,nat1,index_all,x1,mass);
    reset_x(is_lsq,index_lsq,nat1,index_all,x2,mass);
    do_fit(nat1,mass,x1,x2);
  }
  
  bRMS = opt2bSet("-or",NFILE,fnm);
  if (bRMS) {
    fp = xvgropen(opt2fn("-or",NFILE,fnm),"RMSD","Conf","(nm)",oenv);
    xvgr_legend(fp,asize(leg),leg,oenv);
    printf("Select group for RMSD calculation:\n");
    get_index(&atoms,opt2fn_null("-n",NFILE,fnm),1,&isize,&index,&grpname);
    printf("You selected group %s, containing %d atoms\n",grpname,isize);
    rms1 = rmsdev_ind(isize,index,mass,x1,x2);  
    fprintf(stderr,"RMSD between input conformations is %g nm\n",rms1);
  }
  
  snew(dummy,nat1);
  for(i=0; (i<nat1); i++)
    dummy[i] = i;
  status = open_trx(ftp2fn(efTRX,NFILE,fnm),"w");
  
  for(i=0; (i<ninterm); i++) {
    fac = dointerp(nat1,x1,x2,xx,i,ninterm,first,last);
    write_trx(status,nat1,dummy,&atoms,i,fac,box,xx,NULL,NULL);
    if (bRMS) {
      rms1 = rmsdev_ind(isize,index,mass,x1,xx);
      rms2 = rmsdev_ind(isize,index,mass,x2,xx);
      fprintf(fp,"%10g  %10g  %10g\n",fac,rms1,rms2);
    }
  }
  
  close_trx(status); 
  
  if (bRMS) {
    ffclose(fp);
    do_view(oenv,opt2fn("-or",NFILE,fnm),"-nxy");
  }
  
  thanx(stderr);

  return 0;
}
