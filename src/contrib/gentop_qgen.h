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

#ifndef _gentop_qgen_h
#define _gentop_qgen_h
	
#include <stdio.h>
#include "grompp.h"
	
enum { eQGEN_OK, eQGEN_NOTCONVERGED, eQGEN_ERROR, eQGEN_NR };

extern int generate_charges_sm(FILE *fp,void *eem,
			       t_atoms *atoms,rvec x[],
			       real tol,int maxiter,void *atomprop,
			       real qtotref,int eemtype,real hfac,
			       int slater_max);

extern int generate_charges(FILE *fp,char *molname,
			    int eemtype,t_atoms *atoms,rvec x[],
			    real tol,int maxiter,
			    void *atomprop,real qtotref,real hfac);

extern void qgen_message(FILE *fp,int eQGEN);

#endif
