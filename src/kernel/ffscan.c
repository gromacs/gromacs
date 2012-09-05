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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
#include "index.h"
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
    "performed if necessary. Also, like in [TT]mdrun[tt] table functions can be used",
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
    { efDAT, "-ga",     "genalg",   ffOPTRD },
    { efGRO, "-c",      "junk",     ffWRITE },
    { efEDR, "-e",      "junk",     ffWRITE },
    { efTRN, "-o",      "junk",     ffWRITE }
  };
#define NFILE asize(fnm)

  /* Command line options !                         */
  static t_ffscan ff = {
    /* tol      */   0.1,
    /* f_max    */ 100.0,
    /* npow     */  12.0,
    /* epot     */   0.0,
    /* fac_epot */   1.0,
    /* fac_pres */   0.1,
    /* fac_msf  */   0.1,
    /* pres     */   1.0,
    /* molsize  */   1,
    /* nmol     */   1,
    /* bComb    */   TRUE,
    /* bVerbose */   FALSE,
    /* bLogEps  */   FALSE
  };
  static char *loadx=NULL,*loady=NULL,*loadz=NULL;
  static t_pargs pa[] = {
    { "-tol",   FALSE, etREAL, {&ff.tol},   "Energy tolerance (kJ/mol) (zero means everything is printed)" },
    { "-fmax",  FALSE, etREAL, {&ff.f_max},  "Force tolerance (zero means everything is printed)" },
    { "-comb",  FALSE, etBOOL, {&ff.bComb},    "Use combination rules" },
    { "-npow",  FALSE, etREAL, {&ff.npow},     "Power for LJ in case of table use" },
    { "-logeps",FALSE, etBOOL, {&ff.bLogEps},  "Use a logarithmic scale for epsilon" },
    { "-v",     FALSE, etBOOL, {&ff.bVerbose}, "Be loud and noisy" },
    { "-epot",  FALSE, etREAL, {&ff.epot},     "Target energy (kJ/mol)" },
    { "-fepot", FALSE, etREAL, {&ff.fac_epot}, "Factor for scaling energy violations (0 turns energy contribution off)" },
    { "-pres",  FALSE, etREAL, {&ff.pres},     "Value for reference pressure" },
    { "-fpres", FALSE, etREAL, {&ff.fac_pres}, "Factor for scaling pressure violations (0 turns pressure contribution off)" },
    { "-fmsf",  FALSE, etREAL, {&ff.fac_msf},  "Factor for scaling mean square force violations (0 turns MSF contribution off)" },
    { "-molsize",FALSE,etINT,  {&ff.molsize},  "Number of atoms per molecule" },
    { "-nmol",  FALSE, etINT,  {&ff.nmol},     "Number of molecules (Epot is divided by this value!)" }
  };
#define NPA asize(pa)
  unsigned  long Flags = 0;
  gmx_edsam_t ed=NULL;
  FILE      *fplog;

  ivec ddxyz = { 1,1,1 };

  cr = init_par(&argc,&argv);
  
  ff.bVerbose = ff.bVerbose && MASTER(cr);
#if 0
  snew(ed,1);
  ed->eEDtype=eEDnone;
#endif
  
  if (MASTER(cr))
    CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_BE_NICE,NFILE,fnm,
		    NPA,pa,asize(desc),desc,0,NULL);

  if (ff.npow <= 6.0)
    gmx_fatal(FARGS,"Can not have repulsion with smaller exponent than 6");
  if (ff.nmol < 1)
    gmx_fatal(FARGS,"Can not fit %d molecules",ff.nmol);
    
  gmx_log_open(ftp2fn(efLOG,NFILE,fnm),cr,FALSE,0,&fplog);

  if (MASTER(cr)) {
    CopyRight(fplog,argv[0]);
    please_cite(fplog,"Lindahl2001a");
    please_cite(fplog,"Berendsen95a");
  }
  
  set_ffvars(&ff);
  
  Flags = (Flags | MD_FFSCAN);

  mdrunner(fplog,cr,NFILE,fnm,ff.bVerbose,FALSE,
	   ddxyz,0,0,0,loadx,loady,loadz,1,
	   0,0,Flags);
  if (gmx_parallel_env_initialized())
    gmx_finalize(cr);

  gmx_log_close(fplog);

  return 0;
}

