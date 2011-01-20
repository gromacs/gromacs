/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "statutil.h"
#include "viewit.h"
#include "string2.h"
#include "filenm.h"
#include "macros.h"
#include "gmx_fatal.h"

static const int can_view_ftp[] = { 0,
  efEPS,           efXPM,         efXVG,          efPDB };
#define NVIEW asize(can_view_ftp)
static const char* view_program[] = { NULL,
  "ghostview",    "display",      NULL,           "xterm -e rasmol" };

int can_view(int ftp)
{
  int i;
  
  for(i=1; i<NVIEW; i++)
    if ( ftp == can_view_ftp[i] )
      return i;
  
  return 0;
}

void do_view(const output_env_t oenv,const char *fn, const char *opts)
{
  char buf[STRLEN], env[STRLEN];
  const char *cmd;
  int ftp, n;
  
  if (output_env_get_view(oenv) && fn) {
    if (getenv("DISPLAY") == NULL) {
      fprintf(stderr,"Can not view %s, no DISPLAY environment variable.\n",fn);
    } else {
      ftp=fn2ftp(fn);
      sprintf(env, "GMX_VIEW_%s", ftp2ext(ftp));
      upstring(env);
      switch(ftp) {
      case efXVG:
	if ( ! (cmd=getenv(env)) ) {
	  if ( getenv("XMGR") )
	    cmd="xmgr";
	  else
	    cmd="xmgrace";
	}
	break;
      default:
      if ( (n=can_view(ftp)) ) {
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
#ifdef GMX_NO_SYSTEM
        printf("Warning-- No calls to system(3) supported on this platform.");
        printf("Warning-- Skipping execution of 'system(\"%s\")'.", buf);
#else
	if( 0 != system(buf) )
	{
	  gmx_fatal(FARGS,"Failed executing command: %s",buf);
	}
#endif
      }
    }
  }
}

void view_all(const output_env_t oenv,int nf, t_filenm fnm[])
{
  int i;
  
  for(i=0; i<nf; i++)
    if ( can_view(fnm[i].ftp) && is_output(&(fnm[i])) && 
	 ( ! is_optional(&(fnm[i])) || is_set(&(fnm[i])) ) )
      do_view(oenv,fnm[i].fns[0], NULL);
}
