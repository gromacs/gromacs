/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
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
 * Great Red Owns Many ACres of Sand 
 */

#ifndef _do_fit_h
#define _do_fit_h

static char *SRCID_do_fit_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) do_fit.h 1.8 2/2/97"
#endif /* HAVE_IDENT */
extern void do_fit(int natoms,real *w_rls,rvec *xp,rvec *x);
/* Do a least squares fit of x to xp. Atoms which have zero mass
 * (w_rls[i]) are not take into account in fitting.
 * This makes is possible to fit eg. on Calpha atoms and orient
 * all atoms. The routine only fits the rotational part,
 * therefore both xp and x should be centered round the origin.
 */

extern void reset_x(int ncm,atom_id ind_cm[],
		    int nrms,atom_id ind_rms[],rvec x[],real mass[]);
/* Put the center of mass of atoms in the origin 
 * The center of mass is computed from the index ind_cm, while
 * the atoms in ind_rms are reset.
 */

#endif	/* _do_fit_h */
