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
#ifdef XMDRUN
    "xmdrun is the experimental MD program. New features are tested in this",
    "program before being implemented in the default mdrun. Currently under",
    "investigation are: polarizibility, glass simulations, ",
    "Free energy perturbation, X-Ray bombardments",
    "and parallel independent simulations."
#else
    "The mdrun program performs Molecular Dynamics simulations.",
#endif
    "It reads the run input file ([TT]-s[tt]) and distributes the",
    "topology over nodes if needed. The coordinates are passed",
    "around, so that computations can begin.",
    "First a neighborlist is made, then the forces are computed.",
    "The forces are globally summed, and the velocities and",
    "positions are updated. If necessary shake is performed to constrain",
    "bond lengths and/or bond angles.",
    "Temperature and Pressure can be controlled using weak coupling to a",
    "bath.[PAR]",
    "mdrun produces at least three output file, plus one log file",
    "([TT]-g[tt]) per node.",
    "The trajectory file ([TT]-o[tt]), contains coordinates, velocities and",
    "optionally forces.",
    "The structure file ([TT]-c[tt]) contains the coordinates and",
    "velocities of the last step.",
    "The energy file ([TT]-e[tt]) contains energies, the temperature,",
    "pressure, etc, a lot of these things are also printed in the log file",
    "of node 0.",
    "Optionally coordinates can be written to a compressed trajectory file",
    "([TT]-x[tt]).[PAR]",
    "When running in parallel with PVM or an old version of MPI the",
    "[TT]-np[tt] option must be given to indicate the number of",
    "nodes.[PAR]",
    "The option [TT]-dgdl[tt] is only used when free energy perturbation is",
    "turned on.[PAR]",
    "With [TT]-rerun[tt] an input trajectory can be given for which ",
    "forces and energies will be (re)calculated. Neighbor searching will be",
    "performed for every frame, unless [TT]nstlist[tt] is zero",
    "(see the [TT].mdp[tt] file).[PAR]",
    "ED (essential dynamics) sampling is switched on by using the [TT]-ei[tt]",
    "flag followed by an [TT].edi[tt] file.",
    "The [TT].edi[tt] file can be produced using options in the essdyn",
    "menu of the WHAT IF program. mdrun produces a [TT].edo[tt] file that",
    "contains projections of positions, velocities and forces onto selected",
    "eigenvectors.[PAR]",
    "The options [TT]-pi[tt], [TT]-po[tt], [TT]-pd[tt], [TT]-pn[tt] are used",
    "for potential of mean force calculations and umbrella sampling.",
    "See manual.[PAR]",
    "When mdrun receives a TERM signal, it will set nsteps to the current",
    "step plus one. When mdrun receives a USR1 signal, it will set nsteps",
    "to the next multiple of nstxout after the current step.",
    "In both cases all the usual output will be written to file.",
    "When running with MPI, a signal to one of the mdrun processes",
    "is sufficient, this signal should not be sent to mpirun or",
    "the mdrun process that is the parent of the others."
  };
  t_commrec    *cr;
  static t_filenm fnm[] = {
    { efTPX, NULL,      NULL,       ffREAD },
    { efTRN, "-o",      NULL,       ffWRITE },
    { efXTC, "-x",      NULL,       ffOPTWR },
    { efSTO, "-c",      "confout",  ffWRITE },
    { efENX, "-e",      "ener",     ffWRITE },
    { efLOG, "-g",      "md",       ffWRITE },
    { efXVG, "-dgdl",   "dgdl",     ffOPTWR },
    { efTRX, "-rerun",  "rerun",    ffOPTRD },
    { efEDI, "-ei",     "sam",      ffOPTRD },
    { efEDO, "-eo",     "sam",      ffOPTWR },
#ifdef XMDRUN
    { efGCT, "-j",      "wham",     ffOPTRD },
    { efGCT, "-jo",     "bam",      ffOPTRD },
    { efXVG, "-ffout",  "gct",      ffOPTWR },
    { efXVG, "-devout", "deviatie", ffOPTWR },
    { efXVG, "-runav",  "runaver",  ffOPTWR },
#endif
    { efPPA, "-pi",     "pull",     ffOPTRD },
    { efPPA, "-po",     "pullout",  ffOPTWR },
    { efPDO, "-pd",     "pull",     ffOPTWR },
    { efNDX, "-pn",     "pull",     ffOPTRD },
  };
#define NFILE asize(fnm)

  /* Command line options ! */
  static bool bVerbose     = FALSE;
  static bool bCompact     = TRUE;
#ifdef XMDRUN
  static bool bMultiSim    = FALSE;
  static bool bGlas        = FALSE;
  static bool bIonize      = FALSE;
#endif
#ifdef PARALLEL
  static int  nnodes=1;
#endif
  static int  nDLB=0,nstepout=10;
  static t_pargs pa[] = {
#ifdef PARALLEL
    { "-np",      FALSE, etINT, {&nnodes},
      "Number of nodes, must be the same as used for grompp" },
#endif
    { "-v",       FALSE, etBOOL,{&bVerbose}, "Be loud and noisy" },
    { "-compact", FALSE, etBOOL,{&bCompact}, "Write a compact log file" },
#ifdef XMDRUN
    { "-multi",   FALSE, etBOOL,{&bMultiSim}, "Do multiple simulations in parallel (only with -np > 1)" },
    { "-glas",    FALSE, etBOOL,{&bGlas},
      "Do glass simulation with special long range corrections" },
    { "-ionize",  FALSE, etBOOL,{&bIonize},
      "Do a simulation including the effect of an X-Ray bombardment on your system" },
#endif
    { "-dlb",     FALSE, etINT, {&nDLB},
      "HIDDENUse dynamic load balancing every ... step. BUGGY do not use" },
    { "-stepout", FALSE, etINT, {&nstepout},
      "HIDDENFrequency of writing the remaining runtime" }
  };
  t_edsamyn edyn;
  unsigned long Flags;
  
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

#ifdef XMDRUN
  if (bMultiSim && PAR(cr))
    cr = init_msim(cr,NFILE,fnm);
#endif

  if (MASTER(cr)) {
    CopyRight(stdlog,argv[0]);
    please_cite(stdlog,"Berendsen95a");
  }
  
  if (opt2bSet("-ei",NFILE,fnm)) 
    ed_open(NFILE,fnm,&edyn);
    
  Flags = opt2bSet("-rerun",NFILE,fnm) ? MD_RERUN : 0;
#ifdef XMDRUN
  Flags = (Flags | 
	   (bIonize   ? MD_IONIZE   : 0) |
	   (bMultiSim ? MD_MULTISIM : 0) |
	   (bGlas     ? MD_GLAS     : 0));
#endif

  mdrunner(cr,NFILE,fnm,bVerbose,bCompact,nDLB,FALSE,nstepout,&edyn,Flags);

#ifdef USE_MPI
  if (gmx_parallel)
    MPI_Finalize();
#endif  

  return 0;
}

