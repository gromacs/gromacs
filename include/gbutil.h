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
static char *SRCID_gbutil_h = "$Id$";

extern void rotate_conf(int natom,rvec *x,rvec *v,real alfa, real beta,real gamma);
/*rotate() rotates a configuration alfa degrees around the x_axis and beta degrees around the y_axis, *v can be NULL */

extern void orient(int natom,rvec *x,rvec *v, rvec angle,matrix box);
/*orient() rotates a configuration until the largest atom-atom distance is 
 *placed along the z-axis and the second largest distance is placed along
 *the y-axis. Finally the third longest distance is placed along the x-axis
 */

extern void genconf(t_atoms *atoms,rvec *x,real *r,matrix box,ivec n_box);
/*genconf() generates a new configuration by adding boxes*/
extern void gen_box(int NTB,int natoms,rvec *x, matrix box,rvec box_space,
		    bool bCenter);
/* gen_box() generates a box around a configuration, box_space is optional 
 * extra space around it. If NTB = 1 then a truncated octahedon will be 
 * generated (don't!) if bCenter then coordinates will be centered in the 
 * genereated box
 */
