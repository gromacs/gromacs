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
 * S  C  A  M  O  R  G
 */

#include "typedefs.h"
#include "names.h"

char *eblock_names[ebNR+1]=
{
  "CGS","MOLS","SBLOCKS",NULL
};

char *eboxtype_names[ebtNR+1]=
{
  "Cubic","Rectangular", NULL
};

char *ens_names[enNR+1]=
{
  "Grid","Simple", NULL
};

char *ei_names[eiNR+1]=
{ 
  "md", "steep", "cg", "ld", NULL 
};

char *bool_names[BOOL_NR+1]=
{
  "FALSE","TRUE", NULL
};

char *yesno_names[BOOL_NR+1]=
{
  "no","yes", NULL
};

char *ptype_str[eptNR+1] = {
  "Atom", "Nucleus", "Shell", "Bond", "Dummy", NULL
};

char *eel_names[eelNR+1] = {
  "Twin-Range", "Reaction-Field", "Generalized-Reaction-Field",
  "PME", "PPPM", "Switch", "Shift", NULL
};

char *eshake_names[estNR+1] = {
  "Lincs", "Shake", NULL
};

char *egrp_nm[egNR+1] = { 
  "Coul-SR","LJ","Buck", "Coul-LR", "Coul-14", "LJ-14", NULL
};

