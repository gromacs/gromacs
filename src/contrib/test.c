/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.99_development_20071104
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "gromacs/commandline/pargs.h"
#include "copyrite.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/fileio/xvgr.h"
#include "viewit.h"
#include "gromacs/fileio/pdbio.h"
#include "macros.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/units.h"
#include "names.h"
#include "txtdump.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/fileio/confio.h"

real pot(real x,real qq,real c6,real c12)
{
  return c12*pow(x,-12)-c6*pow(x,-6)+qq*ONE_4PI_EPS0/x;
}

real dpot(real x,real qq,real c6,real c12)
{
  return -(12*c12*pow(x,-13)-6*c6*pow(x,-7)+qq*ONE_4PI_EPS0/sqr(x));
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "Plot the potential"
  };
  static real c6=1.0e-3,c12=1.0e-6,qi=1,qj=2,sig=0.3,eps=1,sigfac=0.7;
  t_pargs pa[] = {
    { "-c6",   FALSE,  etREAL,  {&c6},  "c6"   },
    { "-c12",  FALSE,  etREAL,  {&c12}, "c12"  },
    { "-sig",  FALSE,  etREAL,  {&sig}, "sig"  },
    { "-eps",  FALSE,  etREAL,  {&eps}, "eps"  },
    { "-qi",   FALSE,  etREAL,  {&qi},  "qi"   },
    { "-qj",   FALSE,  etREAL,  {&qj},  "qj"   },
    { "-sigfac", FALSE, etREAL, {&sigfac}, "Factor in front of sigma for starting the plot" }
  };
  t_filenm fnm[] = {
    { efXVG, "-o", "potje", ffWRITE }
  };
#define NFILE asize(fnm)

  FILE      *fp;
  int       i;
  real      qq,x,oldx,minimum,mval,dp[2],pp[2];
  int       cur=0;
#define next (1-cur)
  
  /* CopyRight(stdout,argv[0]);*/
  parse_common_args(&argc,argv,PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),
		    desc,0,NULL);

  if (opt2parg_bSet("-sig",asize(pa),pa) ||
      opt2parg_bSet("-eps",asize(pa),pa)) {
    c6  = 4*eps*pow(sig,6);
    c12 = 4*eps*pow(sig,12);
  }
  else if ((c6 != 0) && (c12 != 0)) {
    sig = pow(c12/c6,1.0/6.0);
    eps = c6*c6/(4*c12);
  }
  else {
    sig = eps = 0;
  }
  printf("c6    = %12.5e, c12     = %12.5e\n",c6,c12);
  printf("sigma = %12.5f, epsilon = %12.5f\n",sig,eps);
  qq = qi*qj;
      
  fp = xvgropen(ftp2fn(efXVG,NFILE,fnm),"Potential","r (nm)","E (kJ/mol)");
  if (sig == 0)
    sig=0.25;
  minimum = -1;
  mval    = 0;
  oldx    = 0;
  for(i=0; (i<100); i++) {
    x    = sigfac*sig+sig*i*0.02;
    dp[next] = dpot(x,qq,c6,c12);
    fprintf(fp,"%10g  %10g  %10g\n",x,pot(x,qq,c6,c12),
	    dp[next]);
    if ((i > 0) && (dp[cur]*dp[next] < 0)) {
      minimum = oldx + dp[cur]*(x-oldx)/(dp[cur]-dp[next]);
      mval    = pot(minimum,qq,c6,c12);
      /*fprintf(stdout,"dp[cur] = %g, dp[next] = %g  oldx = %g, dx = %g\n",
	dp[cur],dp[next],oldx,x-oldx);*/
      printf("Minimum at r = %g (nm). Value = %g (kJ/mol)\n",
	      minimum,mval);
    }
    cur = next;
    oldx = x;
      
  }
  gmx_ffclose(fp);
  
  do_view(ftp2fn(efXVG,NFILE,fnm),NULL);

  gmx_thanx(stderr);
	       
  return 0;
}


