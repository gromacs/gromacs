/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Glycine aRginine prOline Methionine Alanine Cystine Serine
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
