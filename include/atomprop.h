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

static char *SRCID_atomprop_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

extern real get_mass(char *resnm, char *atomnm);
/* search the mass belonging to residue and atom,
   note that the longest match is returned */

extern real get_vdw(char *resnm, char *atomnm, real default_r);
/* search the vdw radius belonging to residue and atom ,
   note that the longest match is returned */
