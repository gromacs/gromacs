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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_viewit_c = "$Id$";

#include "statutil.h"
#include "viewit.h"

void do_view(char *fn, char *opts)
{
  char buf[256], *cmd;

  if (bDoView() && fn) {
    if (getenv("DISPLAY") == NULL) {
      fprintf(stderr,"Can not view %s, no DISPLAY environment variable.\n",
	      fn);
    }
    else {
      switch(fn2ftp(fn)) {
      case efEPS:
	if ( ! (cmd=getenv("GMX_VIEW_EPS")) )
	  cmd="ghostview";
	break;
      case efXPM:
	if ( ! (cmd=getenv("GMX_VIEW_XPM")) )
	  cmd="xv";
	break;
      case efXVG:
	if ( ! (cmd=getenv("GMX_VIEW_XVG")) ) {
	  if ( getenv("XMGRACE") )
	    cmd="xmgrace";
	  else
	    cmd="xmgr";
	}
	break;
      case efPDB:
	if ( ! (cmd=getenv("GMX_VIEW_PDB")) )
	  cmd="xterm -e rasmol";
	break;
      default:
	fprintf(stderr,"Don't know how to view file %s",fn);
	return;
      }
      if ( strlen(cmd) ) {
	sprintf(buf,"%s %s %s &",cmd,fn,opts ? opts : "");
	fprintf(stderr,"executing '%s'\n",buf);
	system(buf);
      }
    }
  }
}

