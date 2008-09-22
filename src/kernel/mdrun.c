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
#include "checkpoint.h"

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
    "The user only needs to change this when the system is very inhomogeneous.",
    "Dynamic load balancing is set with the option [TT]-dlb[tt],",
    "which can give a significant performance improvement,",
    "especially for inhomogeneous systems. The only disadvantage of",
    "dynamic load balancing is that runs are no longer binary reproducible,",
    "but in most cases this is not important.",
    "By default the dynamic load balancing is automatically turned on",
    "when the measured performance loss due to load imbalance is 5% or more.",
    "At low parallelization these are the only important options",
    "for domain decomposition.",
    "At high parallelization the options in the next two sections",
    "could be important for increasing the performace.",
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
    "[BR]",
    "Option [TT]-rdd[tt] can be used to set the required maximum distance",
    "for inter charge-group bonded interactions.",
    "Communication for two-body bonded interactions below the non-bonded",
    "cut-off distance always comes for free with the non-bonded communication.",
    "Atoms beyond the non-bonded cut-off are only communicated when they have",
    "missing bonded interactions; this means that the extra cost is minor",
    "and nearly indepedent of the value of [TT]-rdd[tt].",
    "With dynamic load balancing option [TT]-rdd[tt] also sets",
    "the lower limit for the domain decomposition cell sizes.",
    "By default [TT]-rdd[tt] is determined by mdrun based on",
    "the initial coordinates. The chosen value will be a balance",
    "between interaction range and communication cost.",
    "[BR]",
    "When inter charge-group bonded interactions are beyond",
    "the bonded cut-off distance, mdrun terminates with an error message.",
    "For pair interactions and tabulated bonds",
    "that do not generate exclusions, this check can be turned off",
    "with the option [TT]-noddcheck[TT].",
    "[BR]",
    "When constraints are present, option [TT]-rcon[tt] influences",
    "the cell size limit as well.",
    "Atoms connected by NC constraints, where NC is the LINCS order plus 1,",
    "should not be beyond the smallest cell size. A error message is",
    "generated when this happens and the user should change the decomposition",
    "or decrease the LINCS order and increase the number of LINCS iterations.",
    "By default mdrun estimates the minimum cell size required for P-LINCS",
    "in a conservative fashion. For high parallelization it can be useful",
    "to set the distance required for P-LINCS with the option [TT]-rcon[tt].",
    "[BR]",
    "The [TT]-dds[tt] option sets the minimum allowed x, y and/or z scaling",
    "of the cells with dynamic load balancing. mdrun will ensure that",
    "the cells can scale down by at least this factor. This option is used",
    "for the automated spatial decomposition (when not using [TT]-dd[tt])",
    "as well as for determining the number of grid pulses, which in turn",
    "sets the minimum allowed cell size. Under certain circumstances",
    "the value of [TT]-dds[tt] might need to be adjusted to account for",
    "high or low spatial inhomogeneity of the system.",
    "[PAR]",
    "The option [TT]-nosum[tt] can be used to only sum the energies",
    "at every neighbor search step and energy output step.",
    "This can improve performance for highly parallel simulations",
    "where this global communication step becomes the bottleneck.",
    "For a global thermostat and/or barostat the temperature",
    "and/or pressure will also only be updated every nstlist steps.",
    "With this option the energy file will not contain averages and",
    "fluctuations over all integration steps.[PAR]",
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
    "and X-Ray bombardments.",
    "[PAR]",
    "The option [TT]-pforce[tt] is useful when you suspect a simulation",
    "crashes due to too large forces. With this option coordinates and",
    "forces of atoms with a force larger than a certain value will",
    "be printed to stderr.",
    "[PAR]",
    "Checkpoints containing the complete state of the system are written",
    "at regular intervals (option [TT]-cpt[tt]) to the file [TT]-cpo[tt],",
    "unless option [TT]-cpt[tt] is set to -1.",
    "A simulation can be continued by reading the full state from file",
    "with option [TT]-cpi[tt]. This option is intelligent in the way that",
    "if no checkpoint file is found, Gromacs just assumes a normal run and",
	"starts from the first step of the tpr file.",
	"[PAR]",
	"With checkpointing you can also use the option [TT]-append[tt] to",
	"just continue writing to the previous output files. This is not",
	"enabled by default since it is potentially dangerous if you move files,",
	"but if you just leave all your files in place and restart mdrun with",
	"exactly the same command (with options [TT]-cpi[tt] and [TT]-append[tt])",
	"the result will be the same as from a single run. The contents will",
	"be binary identical (unless you use dynamic load balancing),",
	"but for technical reasons there might be some extra energy frames when",
	"using checkpointing (necessary for restarts without appending).",
    "[PAR]",
    "With option [TT]-maxh[tt] a simulation is terminated and a checkpoint",
    "file is written at the first neighbor search step where the run time",
    "exceeds [TT]-maxh[tt]*0.99 hours.",
    "[PAR]",
    "When mdrun receives a TERM signal, it will set nsteps to the current",
    "step plus one. When mdrun receives a USR1 signal, it will stop after",
    "the next neighbor search step (with nstlist=0 at the next step).",
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
    { efCPT, "-cpi",    NULL,       ffOPTRD },
    { efCPT, "-cpo",    NULL,       ffOPTWR },
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
  static bool bDDBondComm  = TRUE;
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
  static char *dddlb_opt[] =
    { NULL, "auto", "no", "yes", NULL };
  static real rdd=0.0,rconstr=0.0,dlb_scale=0.8,pforce=-1;
  static char *ddcsx=NULL,*ddcsy=NULL,*ddcsz=NULL;
  static real cpt_period=15.0,max_hours=-1;
  static bool bAppendFiles=FALSE;
	
  static t_pargs pa[] = {
    { "-pd",      FALSE, etBOOL,{&bPartDec},
      "Use particle decompostion" },
    { "-dd",      FALSE, etRVEC,{&realddxyz},
      "Domain decomposition grid, 0 is optimize" },
    { "-nt",      FALSE, etINT, {&nthreads},
      "HIDDENNumber of threads to start on each node" },
    { "-npme",    FALSE, etINT, {&npme},
      "Number of separate nodes to be used for PME, -1 is guess" },
    { "-ddorder", FALSE, etENUM, {ddno_opt},
      "DD node order" },
    { "-ddcheck", FALSE, etBOOL, {&bDDBondCheck},
      "Check for all bonded interactions with DD" },
    { "-ddbondcomm", FALSE, etBOOL, {&bDDBondComm},
      "HIDDENUse special bonded atom communication when -rdd > cut-off" },
    { "-rdd",     FALSE, etREAL, {&rdd},
      "The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates" },
    { "-rcon",    FALSE, etREAL, {&rconstr},
      "Maximum distance for P-LINCS (nm), 0 is estimate" },
    { "-dlb",     FALSE, etENUM, {dddlb_opt},
      "Dynamic load balancing (with DD)" },
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
    { "-seppot",  FALSE, etBOOL, {&bSepPot},
      "Write separate V and dVdl terms for each interaction type and node to the log file(s)" },
    { "-pforce",  FALSE, etREAL, {&pforce},
      "Print all forces larger than this (kJ/mol nm)" },
    { "-reprod",  FALSE, etBOOL,{&bReproducible},  
      "Try to avoid optimizations that affect binary reproducibility" },
    { "-cpt",     FALSE, etREAL, {&cpt_period},
      "Checkpoint interval (minutes)" },
    { "-append",  FALSE, etBOOL, {&bAppendFiles},
	  "Append to previous output files when restarting from checkpoint" },
    { "-maxh",   FALSE, etREAL, {&max_hours},
      "Terminate after 0.99 times this time (hours)" },
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
  gmx_edsam_t  ed;
  unsigned long Flags, PCA_Flags;
  ivec     ddxyz;
  int      dd_node_order;
  bool     HaveCheckpoint;
  FILE     *fplog,*fptest;
  int      sim_part;
  char     suffix[STRLEN];
	
  cr = init_par(&argc,&argv);

  if (MASTER(cr))
    CopyRight(stderr,argv[0]);

  PCA_Flags = (PCA_KEEP_ARGS | PCA_NOEXIT_ON_ARGS | PCA_CAN_SET_DEFFNM
	       | (MASTER(cr) ? 0 : PCA_QUIET));
  /* Only run niced when not running in parallel */
  if (!gmx_parallel_env)
    PCA_Flags |= PCA_BE_NICE;

  parse_common_args(&argc,argv,PCA_Flags,
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

  /* Check if there is ANY checkpoint file available */	
  sim_part = 1;
  if(opt2bSet("-cpi",NFILE,fnm))
  {
	  sim_part = read_checkpoint_simulation_part(opt2fn("-cpi",NFILE,fnm)) + 1;
	  /* sim_part will now be 1 if no checkpoint file was found */
	  if(sim_part==1)
	  {
		  fprintf(stdout,"No previous checkpoint file present, assuming this is a new run.\n");
	  }
  } 
	
  if (sim_part<=1)
  { 
	  bAppendFiles = FALSE;
  }
	
  if(!bAppendFiles && sim_part > 1)
  {
	  /* This is a continuation run, rename trajectory output files (except checkpoint files) */
	  /* create new part name first (zero-filled) */
	  if(sim_part<10)
		  sprintf(suffix,"part000%d",sim_part);
	  else if(sim_part<100)
		  sprintf(suffix,"part00%d",sim_part);
	  else if(sim_part<1000)
		  sprintf(suffix,"part0%d",sim_part);
	  else
		  sprintf(suffix,"part%d",sim_part);
	  
	  add_suffix_to_output_names(fnm,NFILE,suffix);
	  fprintf(stdout,"Checkpoint file is from part %d, new output files will be suffixed %s.\n",sim_part-1,suffix);
  }	

  Flags = opt2bSet("-rerun",NFILE,fnm) ? MD_RERUN : 0;
  Flags = Flags | (bSepPot       ? MD_SEPPOT       : 0);
  Flags = Flags | (bIonize       ? MD_IONIZE       : 0);
  Flags = Flags | (bGlas         ? MD_GLAS         : 0);
  Flags = Flags | (bPartDec      ? MD_PARTDEC      : 0);
  Flags = Flags | (bDDBondCheck  ? MD_DDBONDCHECK  : 0);
  Flags = Flags | (bDDBondComm   ? MD_DDBONDCOMM   : 0);
  Flags = Flags | (bConfout      ? MD_CONFOUT      : 0);
  Flags = Flags | (!bSumEner     ? MD_NOGSTAT      : 0);
  Flags = Flags | (bReproducible ? MD_REPRODUCIBLE : 0);
  Flags = Flags | (bAppendFiles  ? MD_APPENDFILES  : 0); 
  Flags = Flags | (sim_part>1    ? MD_STARTFROMCPT : 0); 

  
  /* We postpone opening the log file if we are appending, so we can first truncate
   * the old log file and append to the correct position there instead.
   */
  if (MASTER(cr) && !bAppendFiles) 
  {
    fplog = gmx_log_open(ftp2fn(efLOG,NFILE,fnm),cr,!bSepPot,Flags);
    CopyRight(fplog,argv[0]);
    please_cite(fplog,"Hess2008b");
    please_cite(fplog,"Spoel2005a");
    please_cite(fplog,"Lindahl2001a");
    please_cite(fplog,"Berendsen95a");
  }
  else
  {
  	fplog = NULL;
  }
  
  /* Essential dynamics */
  if (opt2bSet("-ei",NFILE,fnm)) {
    /* Open input and output files, allocate space for ED data structure */
    ed = ed_open(NFILE,fnm,cr);
  } else
    ed=NULL;
    
  ddxyz[XX] = (int)(realddxyz[XX] + 0.5);
  ddxyz[YY] = (int)(realddxyz[YY] + 0.5);
  ddxyz[ZZ] = (int)(realddxyz[ZZ] + 0.5);
  
  mdrunner(fplog,cr,NFILE,fnm,bVerbose,bCompact,
	   ddxyz,dd_node_order,rdd,rconstr,
	   dddlb_opt[0],dlb_scale,ddcsx,ddcsy,ddcsz,
	   nstepout,ed,repl_ex_nst,repl_ex_seed,pforce,
	   cpt_period,max_hours,Flags);
  
  if (gmx_parallel_env)
    gmx_finalize(cr);

  if (MULTIMASTER(cr)) {
    thanx(stderr);
  }

	/* Log file has to be closed in mdrunner if we are appending to it (fplog not set here) */
	if (MASTER(cr) && !bAppendFiles) 
	{
		gmx_log_close(fplog);
	}
	
  return 0;
}

