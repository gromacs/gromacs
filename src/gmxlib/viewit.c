/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
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
  "ghostview",    "xv",       "xmgrace",          "xterm -e rasmol" };

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
  char buf[STRLEN], env[20], ext[N_EXT], *cmd, *defopts=NULL;
  int ftp, n;
  
  if (bDoView() && fn) {
    if (getenv("DISPLAY") == NULL) {
      fprintf(stderr,"Can not view %s, no DISPLAY environment variable.\n",fn);
    } else {
      ftp=fn2ftp(fn);
      strncpy(ext, ftp2ext(ftp), N_EXT);
      upstring(ext);
      sprintf(env, "GMX_VIEW_%s", ext);
      if ( (n=can_view(ftp)) ) {
	if ( ! (cmd=getenv(env)) )
	  cmd=view_program[n];
      } else {
	fprintf(stderr,"Don't know how to view file %s",fn);
	return;
      }
      /* Add command line option -nxy for xmgrace */
      if (ftp == efXVG && strcmp(cmd,"xmgrace") == 0)
	defopts = "-nxy";
      if ( strlen(cmd) ) {
	sprintf(buf,"%s %s%s%s%s%s &",
		cmd,
		opts ? opts : "",opts ? " " : "",
		defopts ? defopts : "",defopts ? " " : "",
		fn);
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
