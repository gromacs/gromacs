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

#include <math.h>
#include "typedefs.h"
#include "macros.h"
#include "copyrite.h"
#include "main.h"
#include "statutil.h"
#include "futil.h"
#include "smalloc.h"
#include "string2.h"
#include "edsam.h"
#include "mdrun.h"
#include "strdb.h"
#include "rdgroup.h"
#include "xmdrun.h"
#include "vec.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "The ffscan program performs a single point energy and force calculation",
    "in which the force field is modified. This way a range of parameters can",
    "be changed and tested for reproduction of e.g. quantum chemical or",
    "experimental data. A grid scan over the parameters is done as specified",
    "using command line arguments. All parameters that reproduce the energy",
    "within a given absolute tolerance are printed to a log file.[PAR]",
    "Obviously polarizable models can be used, and shell optimisation is",
    "performed if necessary. Also, like in mdrun table functions can be used",
    "for user defined potential functions.[PAR]",
    "If the option -ga with appropriate file is passed, a genetic algorithm will",
    "be used rather than a grid scan."
  };
  t_commrec    *cr;
  static t_filenm fnm[] = {
    { efTPX, NULL,      NULL,       ffREAD  },
    { efLOG, "-g",      "md",       ffWRITE },
    { efXVG, "-table",  "table",    ffOPTRD },
    { efDAT, "-parm",   "params",   ffREAD  },
    { efDAT, "-ga",     "genalg",   ffOPTRD }
  };
#define NFILE asize(fnm)

  /* Command line options !                         */
  static real tol   = 0.1;
  static real fmax  = 100;
  static real epot  = 0.0;
  static real npow  = 12.0;
  static real ratio = 0.01;
  static bool bComb = TRUE;
  static bool bVerbose = TRUE;
  static bool bLogEps  = FALSE;
  static t_pargs pa[] = {
    { "-tol",   FALSE, etREAL, {&tol},   "Energy tolerance (kJ/mol) (zero means everything is printed)" },
    { "-fmax",  FALSE, etREAL, {&fmax},  "Force tolerance (zero means everything is printed)" },
    { "-epot",  FALSE, etREAL, {&epot},  "Target energy (kJ/mol)" },
    { "-comb",  FALSE, etBOOL, {&bComb}, "Use combination rules" },
    { "-npow",  FALSE, etREAL, {&npow},  "Power for LJ in case of table use" },
    { "-ratio", FALSE, etREAL, {&ratio}, "Ratio for weighing RMS Force and Energy Deviation. Cost is ratio * RMS Force + abs(Energy Deviation). This probably should depend on your system." },
    { "-logeps",FALSE, etBOOL, {&bLogEps},  "Use a logarithmic scale for epsilon" },
    { "-v",     FALSE, etBOOL, {&bVerbose}, "Be loud and noisy" }
  };
  unsigned  long Flags = 0;
  t_edsamyn edyn;
  
  cr = init_par(&argc,&argv);
  
  bVerbose = bVerbose && MASTER(cr);
  edyn.bEdsam=FALSE;
  
  if (MASTER(cr))
    CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_BE_NICE,NFILE,fnm,
		    asize(pa),pa,asize(desc),desc,0,NULL);

  if (npow <= 6.0)
    fatal_error(0,"Can not have repulsion with smaller exponent than 6");

  open_log(ftp2fn(efLOG,NFILE,fnm),cr);

  if (MASTER(cr)) {
    CopyRight(stdlog,argv[0]);
    please_cite(stdlog,"Lindahl2001a");
    please_cite(stdlog,"Berendsen95a");
  }
  
  set_ffvars(tol,epot,npow,bComb,fmax,bLogEps,ratio);
  
  Flags = (Flags | MD_FFSCAN | MD_XMDRUN);

  mdrunner(cr,NFILE,fnm,bVerbose,TRUE,0,1,&edyn,Flags);

  if (gmx_parallel)
    gmx_finalize();

  return 0;
}

