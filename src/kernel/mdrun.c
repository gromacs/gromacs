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
/* afm stuf */
#include "pull.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "The mdrun program performs Molecular Dynamics simulations.",
    "It reads the run input file ([TT]-s[tt]) and distributes the",
    "topology over processors if needed. The coordinates are passed",
    "around, so that computations can begin.",
    "First a neighborlist is made, then the forces are computed.",
    "The forces are globally summed, and the velocities and",
    "positions are updated. If necessary shake is performed to constrain",
    "bond lengths and/or bond angles.",
    "Temperature and Pressure can be controlled using weak coupling to a",
    "bath.[PAR]",
    "mdrun produces at least three output file, plus one log file",
    "([TT]-g[tt]) per processor.",
    "The trajectory file ([TT]-o[tt]), contains coordinates, velocities and",
    "optionally forces.",
    "The structure file ([TT]-c[tt]) contains the coordinates and",
    "velocities of the last step.",
    "The energy file ([TT]-e[tt]) contains energies, the temperature,",
    "pressure, etc, a lot of these things are also printed in the log file",
    "of processor 0.",
    "Optionally coordinates can be written to a compressed trajectory file",
    "([TT]-x[tt]).[PAR]",
    "When running in parallel with PVM or an old version of MPI the",
    "[TT]-np[tt] option must be given to indicate the number of",
    "processors.[PAR]",
    "ED (essential dynamics) sampling is switched on by using the [TT]-ei[tt]",
    "flag followed by an [TT].edi[tt]",
    "file. The [TT].edi[tt] file can be produced using options in the essdyn",
    "menu of the WHAT IF program. mdrun produces a [TT].edo[tt] file that",
    "contains projections of positions, velocities and forces onto selected",
    "eigenvectors.[PAR]",
    "With [TT]-rerun[tt] an input trajectory can be given for which ",
    "forces and energies will be (re)calculated.[PAR]",
    "When mdrun receives a TERM signal it will set nsteps to the current",
    "step plus one, which causes the run to end after one step and write",
    "all the usual output.",
    "When running with MPI, a TERM signal to one of the mdrun processes",
    "is sufficient, this signal should not be sent to mpirun or",
    "the mdrun process that is the parent of the others." 
  };
  t_commrec    *cr;
  static t_filenm fnm[] = {
    { efTPX, NULL, NULL,      ffREAD },
    { efTRN, "-o", NULL,      ffWRITE },
    { efXTC, "-x", NULL,      ffOPTWR },
    { efSTO, "-c", "confout", ffWRITE },
    { efENX, "-e", "ener",    ffWRITE },
    { efLOG, "-g", "md",      ffWRITE },
    { efTRX, "-rerun", "rerun", ffOPTRD },
    /* function "optRerunMDset" (in runner.c) checks if -rerun is specified */
    { efHAT, "-hat","ghat",   ffOPTRD },
    { efEDI, "-ei", "sam",    ffOPTRD },
    { efEDO, "-eo", "sam",    ffOPTWR },
    { efDAT, "-pi","pull",    ffOPTRD },
    { efOUT, "-po","pull",    ffOPTWR },
    { efNDX, "-n", "pull",    ffOPTRD },
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
    { "-np",      FALSE, etINT, &nprocs,
      "Number of processors, must be the same as used for grompp" },
#endif
    { "-v",       FALSE, etBOOL,&bVerbose, "Be loud and noisy" },
    { "-compact", FALSE, etBOOL,&bCompact, "Write a compact log file" },
    { "-dlb",     FALSE, etINT, &nDLB,
      "HIDDENUse dynamic load balancing every ... step. BUGGY do not use" },
    { "-stepout", FALSE, etINT, &nstepout,
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
		    PCA_CAN_SET_DEFFNM | (MASTER(cr) ? 0 : PCA_QUIET),
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
  MPI_Finalize();
#endif  

  return 0;
}

