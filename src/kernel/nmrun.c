/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_nmrun_c = "$Id$";
#include "typedefs.h"
#include "macros.h"
#include "copyrite.h"
#include "main.h"
#include "statutil.h"
#include "futil.h"
#include "mdrun.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "nmrun builds a Hessian matrix from single conformation.",
    "For usual Normal Modes-like calculations, make sure that",
    "the structure provided is properly energy-minimised.",
    "The generated matrix can be diagonalized by g_nmeig."
  };
  t_commrec    *cr;
  t_filenm fnm[] = {
    { efTPX, NULL, NULL,      ffREAD },
    { efMTX, "-m", "hessian", ffWRITE },
    { efLOG, "-g", "nm",      ffWRITE },
  };
#define NFILE asize(fnm)

  /* Command line options ! */
  static bool bVerbose=FALSE,bCompact=TRUE;
  static int  nnodes=1;
  static int  nDLB=0,nstepout=10;
  static t_pargs pa[] = {
    { "-np",      FALSE, etINT, {&nnodes},
      "Number of nodes, must be the same as used for grompp" },
    { "-v",       FALSE, etBOOL,{&bVerbose}, "Verbose mode" },
    { "-compact", FALSE, etBOOL,{&bCompact},
      "Write a compact log file" },
    { "-dlb",     FALSE, etINT, {&nDLB},
      "HIDDENFrequency of dynamic load balancing. BUGGY do not use" },
    { "-stepout", FALSE, etINT, {&nstepout},
      "HIDDENFrequency of writing the remaining runtime" }
  };
  t_edsamyn edyn;
  
  cr = init_par(&argc,&argv);
  bVerbose = bVerbose && MASTER(cr);
  edyn.bEdsam=FALSE;
  
  if (MASTER(cr))
    CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,
		    PCA_KEEP_ARGS | PCA_NOEXIT_ON_ARGS | PCA_BE_NICE | 
		    (MASTER(cr) ? 0 : PCA_QUIET),
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

#ifndef USE_MPI
  if (nnodes>1) 
    fatal_error(0,"GROMACS compiled without MPI support - can't do parallel runs");
#endif
    
  open_log(ftp2fn(efLOG,NFILE,fnm),cr);
  
  if (MASTER(cr)) {
    CopyRight(stdlog,argv[0]);
    please_cite(stdlog,"Berendsen95a");
  }
  
  mdrunner(cr,NFILE,fnm,bVerbose,bCompact,nDLB,TRUE,nstepout,&edyn,0);

#ifdef USE_MPI
  if (gmx_parallel)
    MPI_Finalize();
#endif  

  return 0;
}

