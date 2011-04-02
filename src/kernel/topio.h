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

#ifndef _topio_h
#define _topio_h

#include "typedefs.h"
#include "readir.h"
#include "grompp.h"
#include "gpp_atomtype.h"

extern double check_mol(gmx_mtop_t *mtop,warninp_t wi);
/* Check mass and charge */

extern char **do_top(gmx_bool         bVerbose,
		     const char   *topfile,
		     const char   *topppfile,
		     t_gromppopts *opts,
		     gmx_bool         bZero,
		     t_symtab     *symtab,
		     t_params     plist[],
		     int          *combination_rule,
		     double       *repulsion_power,
		     real         *fudgeQQ,
		     gpp_atomtype_t atype,
		     int          *nrmols,
		     t_molinfo    **molinfo,
		     t_inputrec   *ir,
		     int          *nmolblock,
		     gmx_molblock_t **molblock,
		     gmx_bool         bGB,
		     warninp_t    wi);


/* This routine expects sys->molt[m].ilist to be of size F_NRE and ordered. */
void generate_qmexcl(gmx_mtop_t *sys,t_inputrec *ir);

#endif	/* _topio_h */
