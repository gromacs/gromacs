/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifndef _comp_h
#define _com_h

#ifdef HAVE_IDENT
#ident  "@(#) com.h 1.0 29 Feb 1996"
#endif /* HAVE_IDENT */
#include "typedefs.h"
#include "vec.h"

extern real calc_xcm(rvec x[],        /* coordinates of all atoms        */
		     int ngx,         /* nr of atoms in index group      */
		     atom_id index[], /* the indices of the atoms to use */
		     t_atom atom[],   /* all atoms.                      */
		     rvec xcm);       /* coordinates of center of mass   */
     /* calc_xcm returns the total mass of the molecule */
#endif
