/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
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
#ifdef PARALLEL
  static int  nprocs=1;
#endif
  static int  nDLB=0,nstepout=10;
  static t_pargs pa[] = {
#ifdef PARALLEL
    { "-np",      FALSE, etINT, {&nprocs},
      "Number of processors, must be the same as used for grompp" },
#endif
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
		    PCA_KEEP_ARGS | PCA_NOEXIT_ON_ARGS | PCA_SET_NPRI |
		    (MASTER(cr) ? 0 : PCA_QUIET),
		    TRUE,NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
    
  open_log(ftp2fn(efLOG,NFILE,fnm),cr);
  
  if (MASTER(cr)) {
    CopyRight(stdlog,argv[0]);
    please_cite(stdlog,"Berendsen95a");
  }
  
  mdrunner(cr,NFILE,fnm,bVerbose,bCompact,nDLB,TRUE,nstepout,&edyn);

#ifdef USE_MPI
  if (gmx_parallel)
    MPI_Finalize();
#endif  

  return 0;
}

