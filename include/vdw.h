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
 * Gnomes, ROck Monsters And Chili Sauce
 */
static char *SRCID_vdw_h = "$Id$";

typedef struct {
  char  atomname[10];
  real distance;
} t_vdw;

extern int read_vdw(char *vdwdata,t_vdw **vdw);
/* read_vdw() reads the vdw distances. returns the number of entries read */

extern void write_vdw(char *vdwdata,t_vdw vdw[],int nvdw);
/* writes n vdw distances to a file vdwata */

extern real get_vdw(int nvdw,t_vdw vdw[],char *atom);
/* search the distance belonging to atom */
