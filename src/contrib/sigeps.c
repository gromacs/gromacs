/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Great Red Owns Many ACres of Sand 
 */

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "statutil.h"
#include "copyrite.h"
#include "fatal.h"
#include "xvgr.h"
#include "pdbio.h"
#include "macros.h"
#include "smalloc.h"
#include "vec.h"
#include "pbc.h"
#include "physics.h"
#include "names.h"
#include "txtdump.h"
#include "trnio.h"
#include "symtab.h"
#include "confio.h"

real pot(real x,real qq,real c6,real cn,int npow)
{
  return cn*pow(x,-npow)-c6*pow(x,-6)+qq*ONE_4PI_EPS0/x;
}

real dpot(real x,real qq,real c6,real cn,int npow)
{
  return -(npow*cn*pow(x,-npow-1)-6*c6*pow(x,-7)+qq*ONE_4PI_EPS0/sqr(x));
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "Plot the potential"
  };
  static real c6=1.0e-3,cn=1.0e-6,qi=0,qj=0,sig=0.3,eps=1,sigfac=0.7;
  static int  npow=12;
  t_pargs pa[] = {
    { "-c6",   FALSE,  etREAL,  {&c6},  "c6"   },
    { "-cn",   FALSE,  etREAL,  {&cn},  "constant for repulsion"   },
    { "-pow",  FALSE,  etINT,   {&npow},"power of the repulsion term" },
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
  real      A,B,C;
  int       cur=0;
#define next (1-cur)
  
  /* CopyRight(stdout,argv[0]);*/
  parse_common_args(&argc,argv,PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),
		    desc,0,NULL);

  if (opt2parg_bSet("-sig",asize(pa),pa) ||
      opt2parg_bSet("-eps",asize(pa),pa)) {
    c6  = 4*eps*pow(sig,6);
    cn  = 4*eps*pow(sig,npow);
  }
  else if ((c6 != 0) && (cn != 0)) {
    sig = pow(cn/c6,1.0/(npow-6.0));
    eps = 0.25*c6*pow(sig,-6.0);
  }
  else {
    sig = eps = 0;
  }
  printf("c6    = %12.5e, c%d    = %12.5e\n",c6,npow,cn);
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
    dp[next] = dpot(x,qq,c6,cn,npow);
    fprintf(fp,"%10g  %10g  %10g\n",x,pot(x,qq,c6,cn,npow),
	    dp[next]);
    if (qq != 0) {
      if ((i > 0) && (dp[cur]*dp[next] < 0)) {
	minimum = oldx + dp[cur]*(x-oldx)/(dp[cur]-dp[next]);
	mval    = pot(minimum,qq,c6,cn,npow);
	/*fprintf(stdout,"dp[cur] = %g, dp[next] = %g  oldx = %g, dx = %g\n",
	  dp[cur],dp[next],oldx,x-oldx);*/
	printf("Van der Waals + Coulomb minimum at r = %g (nm). Value = %g (kJ/mol)\n",
	       minimum,mval);
      }
    }
    cur = next;
    oldx = x;
      
  }
  fclose(fp);
  
  if (qq == 0.0) {
    minimum = pow(npow/6.0*pow(sig,npow-6.0),1.0/(npow-6));
    printf("Van der Waals minimum at %g, V = %g\n\n",
	   minimum,pot(minimum,0,c6,cn,npow));
    printf("Fit of Lennard Jones (%d-6) to Buckingham:\n",npow);
    B = npow/minimum;
    C = c6;
    A = 4*eps*pow(sig/minimum,npow)*exp(npow);
    printf("A = %g, B = %g, C = %g\n",A,B,C);
  }
  do_view(ftp2fn(efXVG,NFILE,fnm),NULL);

  thanx(stderr);  
	       
  return 0;
}


