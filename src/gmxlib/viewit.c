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

#include <string.h>
#include "statutil.h"
#include "viewit.h"
#include "string2.h"
#include "filenm.h"
#include "macros.h"

static int can_view_ftp[] = { 0,
  efEPS,           efXPM,         efXVG,          efPDB };
#define NVIEW asize(can_view_ftp)
static char* view_program[] = { NULL,
  "ghostview",    "xv",           NULL,           "xterm -e rasmol" };

int can_view(int ftp)
{
  int i;
  
  for(i=1; i<NVIEW; i++)
    if ( ftp == can_view_ftp[i] )
      return i;
  
  return 0;
}

void do_view(char *fn, char *opts)
{
#define N_EXT 3
  char buf[STRLEN], env[20], ext[N_EXT], *cmd;
  int ftp, n;
  
  if (bDoView() && fn) {
    if (getenv("DISPLAY") == NULL) {
      fprintf(stderr,"Can not view %s, no DISPLAY environment variable.\n",fn);
    } else {
      ftp=fn2ftp(fn);
      strncpy(ext, ftp2ext(ftp), N_EXT);
      upstring(ext);
      sprintf(env, "GMX_VIEW_%s", ext);
      switch(ftp) {
      case efXVG:
	if ( ! (cmd=getenv(env)) ) {
	  if ( getenv("XMGRACE") )
	    cmd="xmgrace";
	  else
	    cmd="xmgr";
	}
	break;
      default:
	if ( n=can_view(ftp) ) {
	  if ( ! (cmd=getenv(env)) )
	    cmd=view_program[n];
	} else {
	  fprintf(stderr,"Don't know how to view file %s",fn);
	  return;
	}
      }
      if ( strlen(cmd) ) {
	sprintf(buf,"%s %s %s &",cmd,opts ? opts : "",fn);
	fprintf(stderr,"Executing '%s'\n",buf);
	system(buf);
      }
    }
  }
}

void view_all(int nf, t_filenm fnm[])
{
  int i;
  
  for(i=0; i<nf; i++)
    if ( can_view(fnm[i].ftp) && is_output(&(fnm[i])) && 
	 ( ! is_optional(&(fnm[i])) || is_set(&(fnm[i])) ) )
      do_view(fnm[i].fn, NULL);
}
