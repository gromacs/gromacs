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

#ifndef _edsam_h
#define _edsam_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "types/edsams.h"

/* this function has to be called before the (LINCS, SHAKE, etc. ) constraints are applied */
extern void prepare_edsam(int step, int start, int homenr, t_commrec *cr, rvec x[],t_edsamyn *edyn);

/* this function implements the edsam, constraints and the .edo monitoring functionality */
extern void do_edsam(FILE *log,t_topology *top,t_inputrec *ir,int step,
		     t_mdatoms *md,int start,int homenr, t_commrec *cr,
                     rvec x[],rvec xold[],rvec f[],matrix box,
                     t_edsamyn *edyn,bool bHave_force);

extern void do_flood(FILE *log, t_commrec *cr, rvec x[],rvec force[], t_edsamyn *edyn, int step);
extern void ed_open(int nfile,t_filenm fnm[],t_edsamyn *edyn, t_commrec *cr);

/* return value is 1 if constraints are switched on, 0 otherwise */
int 
ed_constraints(t_edsamyn *edyn);

extern void init_edsam(FILE *log,t_topology *top,t_inputrec *ir,
		       t_mdatoms *md,int start,int homenr, t_commrec *cr,
		       t_edsamyn *edyn);

extern void do_first_edsam(FILE *log,t_topology *top,
		t_mdatoms *md,int start,int homenr,t_commrec *cr,
			   rvec x[],matrix box, t_edsamyn *edyn,bool bHaveConstr);

extern void finish_edsam(FILE *log,t_topology *top,t_inputrec *ir,
		t_mdatoms *md,int start,int homenr,t_commrec *cr,
			 t_edsamyn *edyn);

#endif	/* _edsam_h */






