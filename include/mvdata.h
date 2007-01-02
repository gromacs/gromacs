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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _mvdata_h
#define _mvdata_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"

extern void ld_data(const t_commrec *cr,int left,int right,
		    t_inputrec *inputrec,
		    t_topology *top,t_state *state);

extern void mv_data(const t_commrec *cr,int left,int right,
		    t_inputrec *inputrec,
		    t_topology *top,t_state *state);

extern void ld_state(const t_commrec *cr,int src,
		     t_state *state);

extern void mv_state(const t_commrec *cr,int dest,
		     t_state *state);

extern void move_cgcm(FILE *log,const t_commrec *cr,rvec cg_cm[]);
		     
extern void move_rvecs(const t_commrec *cr,bool bForward,bool bSum,
		       int left,int right,rvec vecs[],rvec buf[],
		       int shift,t_nrnb *nrnb);

extern void move_x(FILE *log,const t_commrec *cr,
		   int left,int right,rvec x[],t_nrnb *nrnb);
		    
extern void move_f(FILE *log,const t_commrec *cr,
		   int left,int right,rvec f[],rvec fadd[],
		   t_nrnb *nrnb);
		    
#endif	/* _mvdata_h */
