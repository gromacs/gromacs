/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Grunge ROck MAChoS
 */

#ifndef _calch_h
#define _calch_h

static char *SRCID_calch_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) calch.h 1.8 2/2/97"
#endif /* HAVE_IDENT */
#include "typedefs.h"
	
extern void calc_h_pos(int nht, rvec xa[], rvec xh[]);
/*
 *    w.f. van gunsteren, groningen, july 1981 
 *
 *    translated to c d. van der spoel groningen jun 1993
 *    added option 5 jan 95
 *
 *    subroutine genh (nht,nh,na,d,alfa,x)
 *
 *    genh generates cartesian coordinates for hydrogen atoms
 *    using the coordinates of neighbour atoms.
 *
 *    nht      : type of hydrogen attachment (see manual)
 *    xh(1.. ) : atomic positions of the hydrogen atoms that are to be
 *               generated
 *    xa(1..4) : atomic positions of the control atoms i, j and k and l
 *    default bond lengths and angles are defined internally
 */

#endif	/* _calch_h */
