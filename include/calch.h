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
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _calch_h
#define _calch_h

static char *SRCID_calch_h = "$Id$";

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
