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
#include "sysstuff.h"
#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "xvgr.h"
#include "physics.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "nrama.h"
#include "gmx_ana.h"


static void plot_rama(FILE *out,t_xrama *xr)
{
  int i;
  real phi,psi;
  
  for(i=0; (i<xr->npp); i++) {
    phi=xr->dih[xr->pp[i].iphi].ang*RAD2DEG;
    psi=xr->dih[xr->pp[i].ipsi].ang*RAD2DEG;
    fprintf(out,"%g  %g  %s\n",phi,psi,xr->pp[i].label);
  }
}

int gmx_rama(int argc,char *argv[])
{
  const char *desc[] = {
    "[TT]g_rama[tt] selects the [GRK]phi[grk]/[GRK]psi[grk] dihedral combinations from your topology file",
    "and computes these as a function of time.",
    "Using simple Unix tools such as [IT]grep[it] you can select out", 
    "specific residues."
  };

  FILE      *out;
  t_xrama   *xr;
  int       j;
  output_env_t oenv;
  t_filenm  fnm[] = {
    { efTRX, "-f", NULL,  ffREAD },
    { efTPX, NULL, NULL,  ffREAD },
    { efXVG, NULL, "rama",ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,0,NULL,asize(desc),desc,0,NULL,&oenv);

		      
  snew(xr,1);
  init_rama(oenv,ftp2fn(efTRX,NFILE,fnm),ftp2fn(efTPX,NFILE,fnm),xr,3);
  
  out=xvgropen(ftp2fn(efXVG,NFILE,fnm),"Ramachandran Plot","Phi","Psi",oenv);
  xvgr_line_props(out,0,elNone,ecFrank,oenv);
  xvgr_view(out,0.2,0.2,0.8,0.8,oenv);
  xvgr_world(out,-180,-180,180,180,oenv);
  fprintf(out,"@    xaxis  tick on\n@    xaxis  tick major 60\n@    xaxis  tick minor 30\n");
  fprintf(out,"@    yaxis  tick on\n@    yaxis  tick major 60\n@    yaxis  tick minor 30\n");
  fprintf(out,"@ s0 symbol 2\n@ s0 symbol size 0.4\n@ s0 symbol fill 1\n");
  
  j=0;
  do {
    plot_rama(out,xr);
    j++;
  } while (new_data(xr));
  fprintf(stderr,"\n");
  ffclose(out);
  
  do_view(oenv,ftp2fn(efXVG,NFILE,fnm),NULL);
  
  thanx(stderr);
  
  return 0;
}
