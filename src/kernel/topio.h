/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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

#ifndef _topio_h
#define _topio_h

#include "typedefs.h"
#include "readir.h"
#include "grompp.h"

typedef struct {
  int whichmol;
  int nrcopies;
} t_simsystem;

extern real check_mol(t_atoms *atoms);
/* Check mass and charge */

extern void preprocess(char *infile, 
		       char *outfile,
		       char *cpp,
		       char *define,
		       char *include);

extern char **do_top(bool         bVerbose,
		     char         *topfile,
		     char         *topppfile,
		     t_gromppopts *opts,
		     bool         bZero,
		     t_symtab     *symtab,
		     t_params     plist[],
		     int          *combination_rule,
		     real         *repulsion_power,
		     t_atomtype   *atype,
		     int          *nrmols,
		     t_molinfo    **molinfo,
		     t_inputrec   *ir,
		     int          *nsim,
		     t_simsystem  **sims);

void 
generate_qmexcl(t_topology *sys,t_inputrec *ir);

#endif	/* _topio_h */
