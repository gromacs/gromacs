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
 * GROup of MAchos and Cynical Suckers
 */
static char *SRCID_dlg_c = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <xdlghi.h>

int main(int argc, char *argv[])
{
  t_x11 *x11;
  t_dlg *dlg;

  if ((x11=GetX11(&argc,argv))==NULL) {
    fprintf(stderr,"No X!\n");
    exit(1);
  }
  if (argc > 1) {
    dlg=ReadDlg(x11,0,x11->title,x11->fg,x11->bg,argv[1],100,100,TRUE,
		TRUE,NULL,NULL);
    ShowDlg(dlg);
    x11->MainLoop(x11);
  }
  else 
    fprintf(stderr,"Usage: %s [ X options ] infile\n",argv[0]);

  x11->CleanUp(x11);

  return 0;
}
