/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_gbutil_h = "$Id$";

extern int in_box(int NTB,matrix box,rvec x);
/*in_box() returns zero if atom x is inside the box*/

extern void rotate_conf(int natom,rvec *x,rvec *v,real alfa, real beta,real gamma);
/*rotate() rotates a configuration alfa degrees around the x_axis and beta degrees around the y_axis, *v can be NULL */

extern void orient(int natom,rvec *x,rvec *v, rvec angle,matrix box);
/*orient() rotates a configuration until the largest atom-atom distance is 
 *placed along the z-axis and the second largest distance is placed along
 *the y-axis. Finally the third longest distance is placed along the x-axis
 */

extern void genconf(t_atoms *atoms,rvec *x,rvec *v,real *r,matrix box,
		    ivec n_box);
/*genconf() generates a new configuration by adding boxes*/
