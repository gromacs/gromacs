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
 * GRowing Old MAkes el Chrono Sweat
 */

#ifndef	_random_h
#define	_random_h

#ifdef HAVE_IDENT
#ident	"@(#) random.h 1.5 12/16/92"
#endif /* HAVE_IDENT */

#include <typedefs.h>

extern real rando(int *ig);
/* Generate a random number 0 <= r < 1. ig is the (address of) the
 * seed variable.
 */

extern void grp_maxwell(t_block *grp,real tempi[],int nrdf[],int seed,
			t_atoms *atoms,rvec v[]);
/* Generate for each group in grp a temperature. */

extern void maxwell_speed(real tempi,int nrdf,int seed, 
			  t_atoms *atoms, rvec v[]);
/* Generate velocites according to a maxwellian distribution */

extern real calc_cm(FILE *log,int natoms,real mass[],rvec x[],rvec v[],
		    rvec xcm,rvec vcm,rvec acm,matrix L);
/* Calculate the c.o.m. position, velocity, acceleration and the
 * moment of Inertia. Returns the total mass.
 */

extern void stop_cm(FILE *log,int natoms,real mass[],rvec x[],rvec v[]);

#endif	/* _random_h */
