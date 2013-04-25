/*
 * $Id: gentop_qgen.h,v 1.7 2009/01/28 00:04:17 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
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

#ifndef _gentop_qgen_h
#define _gentop_qgen_h

#include <stdio.h>
#include "grompp.h"
#include "poldata.h"
#include "gmx_resp.hpp"
	
enum { eQGEN_OK, eQGEN_NOTCONVERGED, eQGEN_NOSUPPORT, eQGEN_ERROR, eQGEN_NR };

typedef struct gentop_qgen *gentop_qgen_t;

extern gentop_qgen_t 
gentop_qgen_init(gmx_poldata_t pd,t_atoms *atoms,
		 gmx_atomprop_t aps,
		 rvec *x,int eqg_model,real hfac,int qtotal,
		 real epsr);

extern void 
gentop_qgen_done(gentop_qgen_t qgen);

extern int 
generate_charges_sm(FILE *fp,gentop_qgen_t qgen,
		    gmx_poldata_t pd,t_atoms *atoms,rvec x[],
		    real tol,int maxiter,gmx_atomprop_t aps,
		    real hfac,real *chieq);

extern int 
generate_charges(FILE *fp,
		 gentop_qgen_t qgen,
		 gmx_resp_t gr,const char *molname,
		 gmx_poldata_t pd,
		 t_atoms *atoms,rvec x[],
		 real tol,int maxiter,int maxcycle,
		 gmx_atomprop_t aps,real hfac);

extern void 
qgen_message(gentop_qgen_t qgen,int len,char buf[],gmx_resp_t gr);

extern gmx_bool 
bSplitQ(int iModel);

/* The routines below return NOTSET if something is out of the ordinary */
extern int gentop_qgen_get_nzeta(gentop_qgen_t qgen,int atom);

extern int gentop_qgen_get_row(gentop_qgen_t qgen,int atom,int z);

extern double gentop_qgen_get_q(gentop_qgen_t qgen,int atom,int z);

extern double gentop_qgen_get_zeta(gentop_qgen_t qgen,int atom,int z);

#endif
