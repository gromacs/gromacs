/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.99_development_20071104
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "macros.h"
#include "statutil.h"
#include "copyrite.h"

/* 
 * This program is needed to create the files:
 *   options.html
 *   options.tex
 * for the html and latex manuals.
 * It should be ran with the option: -hidden
 */

int main(int argc,char *argv[])
{
  const char *desc[] = {
    "All GROMACS programs have 6 standard options,",
    "of which some are hidden by default:"
  };

  const char *bugs[] = {
    "If the configuration script found Motif or Lesstif on your system, "
    "you can use the graphical interface (if not, you will get an error):[BR]"
    "[TT]-X[tt] gmx_bool [TT]no[tt] Use dialog box GUI to edit command line options",
    
    "When compiled on an SGI-IRIX system, all GROMACS programs have an "
    "additional option:[BR]"
    "[TT]-npri[tt] int [TT]0[tt] Set non blocking priority (try 128)",

    "Optional files are not used unless the option is set, in contrast to "
    "non optional files, where the default file name is used when the "
    "option is not set.",

    "All GROMACS programs will accept file options without a file extension "
    "or filename being specified. In such cases the default filenames will "
    "be used. With multiple input file types, such as generic structure "
    "format, the directory will be searched for files of each type with the "
    "supplied or default name. When no such file is found, or with output "
    "files the first file type will be used.",

    "All GROMACS programs with the exception of [TT]mdrun[tt], "
    "[TT]nmrun[tt] and [TT]eneconv[tt] check if the command line options "
    "are valid.  If this is not the case, the program will be halted.",

    "Enumerated options (enum) should be used with one of the arguments "
    "listed in the option description, the argument may be abbreviated. "
    "The first match to the shortest argument in the list will be selected.",

    "Vector options can be used with 1 or 3 parameters. When only one "
    "parameter is supplied the two others are also set to this value.",

    "For many GROMACS programs, the time options can be supplied in different "
    "time units, depending on the setting of the [TT]-tu[tt] option.",
    
    "All GROMACS programs can read compressed or g-zipped files. There "
    "might be a problem with reading compressed [TT].xtc[tt], "
    "[TT].trr[tt] and [TT].trj[tt] files, but these will not compress "
    "very well anyway.",

    "Most GROMACS programs can process a trajectory with less atoms than "
    "the run input or structure file, but only if the trajectory consists "
    "of the first n atoms of the run input or structure file.",
    
    "Many GROMACS programs will accept the [TT]-tu[tt] option to set the "
    "time units to use in output files (e.g. for [TT]xmgr[tt] graphs or "
    "[TT]xpm[tt] matrices) and in all time options."
  };

  output_env_t oenv=NULL;
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,
		    0,NULL,0,NULL,asize(desc),desc,asize(bugs),bugs,&oenv);
  
  thanx(stderr);
  
  return 0;
}
