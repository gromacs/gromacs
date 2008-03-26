/*
 * $Id$
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "macros.h"
#include "copyrite.h"
#include "main.h"
#include "statutil.h"
#include "smalloc.h"
#include "futil.h"
#include "smalloc.h"
#include "edsam.h"
#include "mdrun.h"
#include "xmdrun.h"
/* afm stuf */
#include "pull.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "The mdrun program is the main computational chemistry engine",
    "within GROMACS. Obviously, it performs Molecular Dynamics simulations,",
    "but it can also perform Stochastic Dynamics, Energy Minimization,",
    "test particle insertion or (re)calculation of energies.",
    "Normal mode analysis is another option. In this case mdrun",
    "builds a Hessian matrix from single conformation.",
    "For usual Normal Modes-like calculations, make sure that",
    "the structure provided is properly energy-minimized.",
    "The generated matrix can be diagonalized by g_nmeig.[PAR]",
    "The mdrun program reads the run input file ([TT]-s[tt])",
    "and distributes the topology over nodes if needed.",
    "mdrun produces at least four output files.",
    "A single log file ([TT]-g[tt]) is written, unless the option",
    "[TT]-seppot[tt] is used, in which case each node writes a log file.",
    "The trajectory file ([TT]-o[tt]), contains coordinates, velocities and",
    "optionally forces.",
    "The structure file ([TT]-c[tt]) contains the coordinates and",
    "velocities of the last step.",
    "The energy file ([TT]-e[tt]) contains energies, the temperature,",
    "pressure, etc, a lot of these things are also printed in the log file.",
    "Optionally coordinates can be written to a compressed trajectory file",
    "([TT]-x[tt]).[PAR]",
    "The option [TT]-dgdl[tt] is only used when free energy perturbation is",
    "turned on.[PAR]",
    "When mdrun is started using MPI with more than 1 node, parallelization",
    "is used. By default domain decomposition is used, unless the [TT]-pd[tt]",
    "option is set, which selects particle decomposition.[PAR]",
    "With domain decomposition, the spatial decomposition can be set",
    "with option [TT]-dd[tt]. By default mdrun selects a good decomposition.",
    "The user only needs to change this when there are errors due to missing",
    "interactions, or when the system is very inhomogeneous.",
    "Dynamic load balancing can be turned on with the option [TT]-dlb[tt],",
    "which can give a significant performance improvement,",
    "especially for inhomogeneous systems. The only disadvantage of",
    "dynamic load balancing is that runs are no longer binary reproducible,",
    "but in most cases this is not important.",
    "For regular atomistic simulations at low or medium parallelization",
    "these are the only important options for domain decomposition.",
    "In other cases, the options in the next two sections",
    "could be essential for allowing a simulation to run or",
    "to optimize performance.",
    "[PAR]",
    "When PME is used with domain decomposition, separate nodes can",
    "be assigned to do only the PME mesh calculation;",
    "this is computationally more efficient starting at about 12 nodes.",
    "The number of PME nodes is set with option [TT]-npme[tt],",
    "this can not be more than half of the nodes.", 
    "By default mdrun makes a guess for the number of PME",
    "nodes when the number of nodes is larger than 11 or performance wise",
    "not compatible with the PME grid x dimension.",
    "But the user should optimize npme. Performance statistics on this issue",
    "are written at the end of the log file.",
    "For good load balancing at high parallelization,",
    "npme should be divisible by the number of PME nodes",
    "[PAR]",
    "This section lists all options that affect the domain decomposition.",
    "All distances required for two-body bonded interactions should be within",
    "the non-bonded cut-off distance; all distances required for multi-body",
    "bonded interactions should be within the smallest cell size.",
    "If this is not the case mdrun terminates with an error message.",
    "For pair interactions and tabulated bonds that do not generate",
    "exclusions, this check can be turned off with the option -noddbc.",
    "By default mdrun sets the minimum cell size to",
    "half the cut-off distance or the maximum distance required",
    "for constraints, whichever is larger.",
    "Option [TT]-rdd[tt] can be used to set the allowed distance for",
    "multi-body bonded interactions or also for two-body bonded interactions",
    "when [TT]-rdd[tt] is larger than the non-bonded cut-off.",
    "With dynamic load balancing option [TT]-rdd[tt] also sets",
    "the lower limit for the domain decomposition cell sizes.",
    "When constraints are present, option [TT]-rcon[tt] influences",
    "the cell size limit as well.",
    "Atoms connected by NC constraints, where NC is the LINCS order plus 1,",
    "should not be beyond the smallest cell size. A error message is",
    "generated when this happens and the user should change the decomposition",
    "or decrease the LINCS order and increase the number of LINCS iterations.",
    "By default mdrun estimates the minimum cell size required for P-LINCS",
    "in a conservative fashion. For high parallelization it can be useful",
    "to set the distance required for P-LINCS with the option [TT]-rcon[tt].",
    "The [TT]-dds[tt] option sets the minimum allowed x, y and/or z scaling",
    "of the cells with dynamic load balancing. mdrun will ensure that",
    "the cells can scale down by at least this factor. This option is used",
    "for the automated spatial decomposition (when not using [TT]-dd[tt])",
    "as well as for determining the number of grid pulses, which in turn",
    "sets the minimum allowed cell size. Under certain circumstances",
    "the value of [TT]-dds[tt] might need to be adjusted to account for",
    "high or low spatial inhomogeneity of the system.",
    "Finally the option [TT]-nosum[tt] can be used to only sum the energies",
    "at every neighbor search step. This can improve performance for highly",
    "parallel simulations, where this global communication step becomes",
    "the bottleneck. This option is currently not allowed when",
    "an extended Hamiltonian thermostat or barostat is used.[PAR]",
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
    "When user-defined potential functions have been selected in the",
    "[TT].mdp[tt] file the [TT]-table[tt] option is used to pass mdrun",
    "a formatted table with potential functions. The file is read from",
    "either the current directory or from the GMXLIB directory.",
    "A number of pre-formatted tables are presented in the GMXLIB dir,",
    "for 6-8, 6-9, 6-10, 6-11, 6-12 Lennard Jones potentials with",
    "normal Coulomb.",
    "When pair interactions are present a separate table for pair interaction",
    "functions is read using the [TT]-tablep[tt] option.[PAR]",
    "When tabulated bonded functions are present in the topology,",
    "interaction functions are read using the [TT]-tableb[tt] option.",
    "For each different tabulated interaction type the table file name is",
    "modified in a different way: before the file extension an underscore is",
    "appended, then a b for bonds, an a for angles or a d for dihedrals",
    "and finally the table number of the interaction type.[PAR]",
    "The options [TT]-pi[tt], [TT]-po[tt], [TT]-pd[tt], [TT]-pn[tt] are used",
    "for potential of mean force calculations and umbrella sampling.",
    "See manual.[PAR]",
    "With [TT]-multi[tt] multiple systems are simulated in parallel.",
    "As many input files are required as the number of systems.",
    "The system number is appended to the run input and each output filename,",
    "for instance topol.tpr becomes topol0.tpr, topol1.tpr etc.",
    "The number of nodes per system is the total number of nodes",
    "divided by the number of systems.",
    "One use of this option is for NMR refinement: when distance",
    "or orientation restraints are present these can be ensemble averaged",
    "over all the systems.[PAR]",
    "With [TT]-replex[tt] replica exchange is attempted every given number",
    "of steps. The number of replicas is set with the [TT]-multi[tt] option,",
    "see above.",
    "All run input files should use a different coupling temperature,",
    "the order of the files is not important. The random seed is set with",
    "[TT]-reseed[tt]. The velocities are scaled and neighbor searching",
    "is performed after every exchange.[PAR]",
    "Finally some experimental algorithms can be tested when the",
    "appropriate options have been given. Currently under",
    "investigation are: polarizability, glass simulations",
    "and X-Ray bombardments.[PAR]",
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
    { efXVG, "-field",  "field",    ffOPTWR },
    { efXVG, "-table",  "table",    ffOPTRD },
    { efXVG, "-tablep", "tablep",   ffOPTRD },
    { efXVG, "-tableb", "table",    ffOPTRD },
    { efTRX, "-rerun",  "rerun",    ffOPTRD },
    { efXVG, "-tpi",    "tpi",      ffOPTWR },
    { efXVG, "-tpid",   "tpidist",  ffOPTWR },
    { efEDI, "-ei",     "sam",      ffOPTRD },
    { efEDO, "-eo",     "sam",      ffOPTWR },
    { efGCT, "-j",      "wham",     ffOPTRD },
    { efGCT, "-jo",     "bam",      ffOPTWR },
    { efXVG, "-ffout",  "gct",      ffOPTWR },
    { efXVG, "-devout", "deviatie", ffOPTWR },
    { efXVG, "-runav",  "runaver",  ffOPTWR },
    { efXVG, "-px",     "pullx",    ffOPTWR },
    { efXVG, "-pf",     "pullf",    ffOPTWR },
    { efMTX, "-mtx",    "nm",       ffOPTWR },
    { efNDX, "-dn",     "dipole",   ffOPTWR }
  };
#define NFILE asize(fnm)

  /* Command line options ! */
  static bool bCart        = FALSE;
  static bool bPPPME       = FALSE;
  static bool bPartDec     = FALSE;
  static bool bDDBondCheck = TRUE;
  static bool bDLB         = FALSE;
  static bool bSumEner     = TRUE;
  static bool bVerbose     = FALSE;
  static bool bCompact     = TRUE;
  static bool bSepPot      = FALSE;
  static bool bGlas        = FALSE;
  static bool bIonize      = FALSE;
  static bool bConfout     = TRUE;
  static bool bReproducible = FALSE;
    
  static int  npme=-1;
  static int  nmultisim=0;
  static int  repl_ex_nst=0;
  static int  repl_ex_seed=-1;
  static int  nstepout=100;
  static int  nthreads=1;

  static rvec realddxyz={0,0,0};
  static char *ddno_opt[ddnoNR+1] =
    { NULL, "interleave", "pp_pme", "cartesian", NULL };
  static real rdd=0.0,rconstr=0.0,dlb_scale=0.8;
  static char *ddcsx=NULL,*ddcsy=NULL,*ddcsz=NULL;

  static t_pargs pa[] = {
    { "-pd",      FALSE, etBOOL,{&bPartDec},
      "Use particle decompostion" },
    { "-dd",      FALSE, etRVEC,{&realddxyz},
      "Domain decomposition grid, 0 is optimize" },
    { "-nt",      FALSE, etINT, {&nthreads},
      "Number of threads to start on each node" },
    { "-npme",    FALSE, etINT, {&npme},
      "Number of separate nodes to be used for PME, -1 is guess" },
    { "-ddorder", FALSE, etENUM, {ddno_opt},
      "DD node order" },
    { "-ddbc",    FALSE, etBOOL, {&bDDBondCheck},
      "Check for all bonded interactions with DD" },
    { "-rdd",     FALSE, etREAL, {&rdd},
      "The minimum distance for DD communication (nm)" },
    { "-rcon",    FALSE, etREAL, {&rconstr},
      "Maximum distance for P-LINCS (nm), 0 is estimate" },
    { "-dlb",     FALSE, etBOOL, {&bDLB},
      "Use dynamic load balancing (only with DD)" },
    { "-dds",     FALSE, etREAL, {&dlb_scale},
      "Minimum allowed dlb scaling of the DD cell size" },
    { "-ddcsx",   FALSE, etSTR, {&ddcsx},
      "HIDDENThe DD cell sizes in x" },
    { "-ddcsy",   FALSE, etSTR, {&ddcsy},
      "HIDDENThe DD cell sizes in y" },
    { "-ddcsz",   FALSE, etSTR, {&ddcsz},
      "HIDDENThe DD cell sizes in z" },
    { "-sum",     FALSE, etBOOL,{&bSumEner},
      "Sum the energies at every step" },
    { "-v",       FALSE, etBOOL,{&bVerbose},  
      "Be loud and noisy" },
    { "-compact", FALSE, etBOOL,{&bCompact},  
      "Write a compact log file" },
    { "-seppot", FALSE, etBOOL,{&bSepPot},
      "Write separate V and dVdl terms for each interaction type and node to the log file(s)" },
    { "-reprod", FALSE, etBOOL,{&bReproducible},  
      "Try to avoid optimizations that affect binary reproducibility" },      
    { "-multi",   FALSE, etINT,{&nmultisim}, 
      "Do multiple simulations in parallel" },
    { "-replex",  FALSE, etINT, {&repl_ex_nst}, 
      "Attempt replica exchange every # steps" },
    { "-reseed",  FALSE, etINT, {&repl_ex_seed}, 
      "Seed for replica exchange, -1 is generate a seed" },
    { "-glas",    FALSE, etBOOL,{&bGlas},
      "Do glass simulation with special long range corrections" },
    { "-ionize",  FALSE, etBOOL,{&bIonize},
      "Do a simulation including the effect of an X-Ray bombardment on your system" },
    { "-confout", FALSE, etBOOL, {&bConfout},
      "HIDDENWrite the last configuration with -c" },
    { "-stepout", FALSE, etINT, {&nstepout},
      "HIDDENFrequency of writing the remaining runtime" }
  };
  t_edsamyn *edyn;
  unsigned long Flags;
  ivec     ddxyz;
  int      dd_node_order;
  FILE     *fplog;

  cr = init_par(&argc,&argv);
  snew(edyn,1);

  if (MASTER(cr))
    CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,
		    PCA_KEEP_ARGS | PCA_NOEXIT_ON_ARGS | PCA_BE_NICE |
		    PCA_CAN_SET_DEFFNM | (MASTER(cr) ? 0 : PCA_QUIET),
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  dd_node_order = nenum(ddno_opt);
  cr->npmenodes = npme;
    
#ifndef GMX_THREADS
  if (nthreads > 1)
    gmx_fatal(FARGS,"GROMACS compiled without threads support - can only use one thread");
#endif

  if (repl_ex_nst != 0 && nmultisim < 2)
    gmx_fatal(FARGS,"Need at least two replicas for replica exchange (option -multi)");

  if (nmultisim > 1 && PAR(cr))
    init_multisystem(cr,nmultisim,NFILE,fnm,TRUE);

  if (MASTER(cr) || bSepPot) {
    fplog = gmx_log_open(ftp2fn(efLOG,NFILE,fnm),cr,!bSepPot);
  } else {
    fplog = NULL;
  }

  if (MASTER(cr)) {
    CopyRight(fplog,argv[0]);
    please_cite(fplog,"Hess2008b");
    please_cite(fplog,"Spoel2005a");
    please_cite(fplog,"Lindahl2001a");
    please_cite(fplog,"Berendsen95a");
  }
  
  if (opt2bSet("-ei",NFILE,fnm)) 
    ed_open(NFILE,fnm,edyn,cr);
    
  Flags = opt2bSet("-rerun",NFILE,fnm) ? MD_RERUN : 0;
  Flags = Flags | (bSepPot       ? MD_SEPPOT       : 0);
  Flags = Flags | (bIonize       ? MD_IONIZE       : 0);
  Flags = Flags | (bGlas         ? MD_GLAS         : 0);
  Flags = Flags | (bPartDec      ? MD_PARTDEC      : 0);
  Flags = Flags | (bDDBondCheck  ? MD_DDBONDCHECK  : 0);
  Flags = Flags | (bDLB          ? MD_DLB          : 0);
  Flags = Flags | (bConfout      ? MD_CONFOUT      : 0);
  Flags = Flags | (!bSumEner     ? MD_NOGSTAT      : 0);
  Flags = Flags | (bReproducible ? MD_REPRODUCIBLE : 0);
    
  ddxyz[XX] = (int)(realddxyz[XX] + 0.5);
  ddxyz[YY] = (int)(realddxyz[YY] + 0.5);
  ddxyz[ZZ] = (int)(realddxyz[ZZ] + 0.5);
  
  mdrunner(fplog,cr,NFILE,fnm,bVerbose,bCompact,
	   ddxyz,dd_node_order,rdd,rconstr,dlb_scale,ddcsx,ddcsy,ddcsz,
	   nstepout,edyn,repl_ex_nst,repl_ex_seed,Flags);
  
  if (gmx_parallel_env)
    gmx_finalize(cr);

  if (MULTIMASTER(cr)) {
    thanx(stderr);
  }
  
  gmx_log_close(fplog);
  
  return 0;
}

