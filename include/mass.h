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
 * GROwing Monsters And Cloning Shrimps
 */
static char *SRCID_mass_h = "$Id$";

typedef struct {
  char atomname[10];
  real mass;
} t_mass;

extern int read_mass(char *massdata,t_mass **mass);
/* read_mass() reads the masses. returns the number of entries read */

extern void write_mass(char *massdata,t_mass mass[],int nmass);
/* writes n masses to a file massdata */

extern real get_mass(int nmass,t_mass mass[],char *atom);
/* search the mass belonging to atom */
