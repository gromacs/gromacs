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
#include "typedefs.h"
#include "xvgr.h"
#include "gstat.h"
#include "copyrite.h"
#include "macros.h"
#include "random.h"
#include "smalloc.h"

int main(int argc,char *argv[])
{
  FILE *fp;
  const char *desc[] = {
    "testac tests the functioning of the GROMACS acf routines"
  };
  static int nframes = 1024;
  static int datatp  = 0;
  static real a=0.02*M_PI;
  output_env_t oenv;
  t_pargs pa[] = {
    { "-np", FALSE, etINT, &nframes,
      "Number of data points" },
    { "-dtp",FALSE, etINT, &datatp,
      "Which data: 0=all 0.0, 1=all 1.0, 2=cos(a t), 3=random, 4=cos(a t)+random, 5=sin(a t)/(a t)" }
  };
  static char *str[] = {
    "all 0.0", 
    "all 1.0",
    "cos(a t)",
    "random", 
    "cos(a t)+random",
    "sin(a t)/(a t)"
  };
  t_filenm fnm[] = {
    { efXVG, "-d", "acf-data", ffWRITE },
    { efXVG, "-c", "acf-corr", ffWRITE },
    { efXVG, "-comb", "acf-comb.xvg", ffWRITE }
  };
#define NFILE asize(fnm)
  int     npargs,i,nlag;
  int     seed=1993;
  real    *data,*data2,x;
  t_pargs *ppa;
  
  CopyRight(stderr,argv[0]);
  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args_r(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
		      NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL,&oenv);
  snew(data,nframes);
  snew(data2,nframes);
  
  fp = xvgropen(opt2fn("-d",NFILE,fnm),"testac","x","y",oenv);
  for(i=0; (i<nframes); i++) {
    x = a*i;
    switch (datatp) {
    case 1:
      data[i] = 1;
      break;
    case 2:
      data[i] = cos(x);
      break;
    case 3:
      data[i] = 2*rando(&seed)-1.0;
      break;
    case 4:
      data[i] = cos(x)+2*rando(&seed)-1.0;
      break;
    case 5:
      if (i==0)
        data[i] = 1;
      else
        data[i] = sin(x)/(x);
      break;
    default:
      /* Data remains 0.0 */
      break;
    }
    fprintf(fp,"%10g  %10g\n",x,data[i]);
    data2[i] = data[i];
  }
  ffclose(fp);
  
  do_autocorr(opt2fn("-c",NFILE,fnm),oenv,str[datatp],
	      nframes,1,&data,a,eacNormal,FALSE);
	      
  nlag = get_acfnout();
  fp = xvgropen(opt2fn("-comb",NFILE,fnm),"testac","x","y",oenv);
  for(i=0; (i<nlag); i++) {
    fprintf(fp,"%10g  %10g  %10g\n",a*i,data2[i],data[i]);
  }
  ffclose(fp);

  do_view(opt2fn("-c",NFILE,fnm),"-nxy");
    
  thanx(stderr);

  return 0;
}
