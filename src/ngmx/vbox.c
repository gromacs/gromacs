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
static char *SRCID_vbox_c = "$Id$";

#include "x11.h"
#include "nmol.h"

t_videobox *init_vbox(t_x11 *x11,t_molwin *mw)
{
  t_videobox *vb;

  snew(vb,1);
  snew(vb->vcr,IDNR-IDBUTNR);
  
  return vb;
}

void done_vbox(t_x11 *x11,t_videobox *vbox)
{
  sfree(vbox->vcr);
  sfree(vbox);
}

