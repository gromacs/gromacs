/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */
static char *SRCID_mdrun_c = "$Id$";

#include "typedefs.h"
#include "macros.h"
#include "copyrite.h"
#include "main.h"
#include "statutil.h"
#include "futil.h"
#include "edsam.h"
#include "mdrun.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "The mdrun program performs Molecular Dynamics simulations.",
    "It reads the run input file (.tpx) file and distributes the",
    "topology over processors if needed. The coordinates are passed",
    "around, so that computations can begin.",
    "First a neighbourlist is made, then the forces are computed.",
    "The forces are globally summed, and the velocities and",
    "positions are updated. If necessary shake is performed to constrain",
    "bond lengths and/or bond angles.",
    "Temperature and Pressure can be controlled using weak coupling to a",
    "bath.[PAR]",
    "A number of environment variables can be set to influence the behaviour",
    "of the mdrun program. Most of these are for debugging purposes, but they",
    "sometimes come in handy when porting the software to an",
    "unsupported platform as well. These environment variables", 
    "are listed elsewhere in the manual.[PAR]",
    "The mdrun produces three output file, plus one log file per processor.",
    "The first file is the trajectory, containing coordinates, velocities",
    "etc. The second file contains the coordinates and velocities at the end",
    "of the run plus the computational box. The third file contains energies.",
    "In the log file from processor 0 the energies, temperature, etc. are printed.[PAR]",
    "When run on a parallel computer or with PVM on a cluster of workstations",
    "the [BB]-np[bb] option must be given to indicate the number of",
    "processors. Note that at current PVM does work, but it may launch",
    "multiple processes on a single processor, which is not very sensible.[PAR]",
    "ED (essential dynamics) sampling is switched on by using the [BB]-ei[bb] flag followed by an .edi",
    "file. The .edi file can be produced using options in the essdyn menu",
    "of the WHAT IF program. mdrun produces a .edo file that contains",
    "projections of positions, velocities and forces onto selected",
    "eigenvectors.[PAR]",
    "With [BB]-rerun[bb] an input trajectory can be given for which ",
    "forces and energies will be (re)calculated.[PAR]",
  };
  char         *lognm=NULL;
  t_commrec    *cr;
  static t_filenm fnm[] = {
    { efTPX, NULL, NULL,      ffREAD },
    { efTRN, "-o", NULL,      ffWRITE },
    { efXTC, "-x", NULL,      ffOPTWR },
    { efGRO, "-c", "confout", ffWRITE },
    { efENX, "-e", "ener",    ffWRITE },
    { efLOG, "-g", "md",      ffWRITE },
    { efTRX, "-rerun", "rerun", ffOPTRD },
    /* function "optRerunMDset" (in runner.c) checks if -rerun is specified */
    { efHAT, "-hat","ghat",   ffOPTRD },
    { efEDI, "-ei", "sam",    ffOPTRD },
    { efEDO, "-eo", "sam",    ffOPTWR }
  };
#define NFILE asize(fnm)

  /* Command line options ! */
  static bool bVerbose=FALSE,bCompact=TRUE;
  static int  nprocs=1,nDLB=0,nstepout=10;
  static t_pargs pa[] = {
    { "-np",      FALSE, etINT, &nprocs,
      "Number of processors, must be the same as used for grompp" },
    { "-v",       FALSE, etBOOL,&bVerbose, "Verbose mode" },
    { "-compact", FALSE, etBOOL,&bCompact,
      "Write a compact log file, i.e. do not write a lot of things which are already in the energy file" },
    { "-dlb",     FALSE, etINT, &nDLB,
      "HIDDENUse dynamic load balancing every ... step. BUGGY do not use" },
    { "-stepout", FALSE, etINT, &nstepout,
      "Frequency of writing the remaining runtime" }
  };
  t_edsamyn edyn;
  
  cr = init_par(argv);
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
  
  if (opt2bSet("-ei",NFILE,fnm)) 
    ed_open(NFILE,fnm,&edyn);
    
  mdrunner(cr,NFILE,fnm,bVerbose,bCompact,nDLB,FALSE,nstepout,&edyn);

#ifdef USE_MPI
  if (nprocs > 1)
    MPI_Finalize();
#endif  

  exit(0);
  
  return 0;
}

