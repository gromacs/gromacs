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
 * GROwing Monsters And Cloning Shrimps
 */
static char *SRCID_plotphi_h = "$Id$";

#ifndef _plot_phi
#define _plot_phi
	
extern void plot_phi(char *fn,rvec box,int natoms,rvec x[],real phi[]);
/* Plot potential (or whatever) in a postscript matrix */

extern void print_phi(char *fn,int natoms,rvec x[],real phi[]);
/* Print to a text file in x y phi format */

extern void plot_qtab(char *fn,int nx,int ny,int nz,real ***qtab);
/* Plot a charge table to a postscript matrix */

#endif
