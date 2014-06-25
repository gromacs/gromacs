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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#include "macros.h"
#include "gromacs/utility/smalloc.h"
#include "xmlio.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"

int cmain(int argc,char *argv[])
{
  int        step,natoms;
  real       t,lambda;
  t_inputrec ir;
  t_topology top;
  matrix     box;
  rvec       *x,*v,*f;
  t_filenm fnm[] = {
    { efTPX, NULL, NULL, ffREAD  },
    { efXML, "-r", NULL, ffREAD  },
    { efXML, "-o", NULL, ffWRITE }
  };  
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,0,NULL,0,NULL,0,NULL);
  
  init_top(&top);
  if (opt2bSet("-r",NFILE,fnm))
    read_xml(opt2fn("-r",NFILE,fnm),&step,&t,&lambda,&ir,
	     box,&natoms,&x,&v,&f,&top);
  else {
    t_tpxheader tpx;
    
    read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&tpx,FALSE);
    snew(x,tpx.natoms);
    snew(v,tpx.natoms);
    f = NULL;
    read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,&ir,
	     box,&natoms,x,v,f,&top);
  }
  /*write_xml(opt2fn("-o",NFILE,fnm),step,t,lambda,&ir,box,natoms,x,v,f,&top);*/
  
  return 0;
}



