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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* This file is completely threadsafe - keep it that way! */


#include <stdarg.h>
#include <ctype.h>
#include "typedefs.h"
#include "string2.h"
#include "gmx_fatal.h"
#include "ffscanf.h"

static int getfld(char **p)
{
  int fld;
  
  fld=0;
  while (isdigit(**p)) fld=(fld*10)+((*((*p)++))-'0');
  return fld;
}

void ffscanf(FILE *in,char *fmt, ...)
{
  va_list ap;
  char *p;
  char buf[STRLEN];
  int i,fld;
  double dval;

  va_start(ap,fmt);
  for (p=fmt; *p; p++) {
    if (*p == '%') {
      p++;
      fld=getfld(&p);
      for(i=0; (i<fld); ) {
	buf[i]=fgetc(in);
	if (buf[i] != '\n') i++;
      }
      buf[fld]='\0';
      switch(*p) {
      case 'd':
	sscanf(buf,"%d",va_arg(ap,int *));
	break;
      case 'f':
	sscanf(buf,"%f",va_arg(ap,float *));
	break;
      case 'F':
	sscanf(buf,"%lf",va_arg(ap,double *));
	break;
      case 'r':
	sscanf(buf,"%lf",&dval);
	*(va_arg(ap,real *)) = dval;
	break;
      default:
	break;
      }
    }
    else
      gmx_fatal(FARGS,"unknown ffscanf format '%c'",*p+1);
  }
  va_end(ap);
}

