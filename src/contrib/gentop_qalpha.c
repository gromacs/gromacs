/*
 * $Id$
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "assert.h"
#include "maths.h"
#include "macros.h"
#include "copyrite.h"
#include "bondf.h"
#include "string2.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "confio.h"
#include "physics.h"
#include "statutil.h"
#include "vec.h"
#include "random.h"
#include "3dview.h"
#include "txtdump.h"
#include "readinp.h"
#include "names.h"
//#include "toppush.h"
//#include "pdb2top.h"
//#include "topexcl.h"
#include "gentop_nm2type.h"
#include "gentop_qalpha.h"

typedef struct {
  char *atom;
  double q,alpha;
} t_qalpha;

typedef struct {
  int nr;
  t_qalpha *qa;
} gentop_qalpha;

gentop_qat rd_q_alpha(char *fn)
/* Read file with charges and polarizabilities */
{
  FILE *fp;
  int    n=0;
  double q,alpha;
  char   name[STRLEN];
  char   line[STRLEN];
  gentop_qalpha *qa=NULL;
  
  snew(qa,1);
  fp = ffopen(fn,"r");
  while (fgets2(line,STRLEN-1,fp) != NULL) {
    strip_comment(line);
    trim(line);
    if (strlen(line) > 0) {
      if (sscanf(line,"%s%lf%lf",name,&q,&alpha) == 3) {
	srenew(qa->qa,n+1);
	qa->qa[n].atom  = strdup(name);
	qa->qa[n].q     = q;
	qa->qa[n].alpha = alpha;
	n++;
      }
      else 
	fprintf(stderr,"Don't understand line '%s'\n",line);
    }
  }
  fclose(fp);
  qa->nr = n;
  
  return (gentop_qat) qa;
}

void dump_q_alpha(FILE *fp,gentop_qat qat)
/* Dump the database for debugging */
{
}

double get_qa_q(char *atom,gentop_qat qat)
/* Return the charge belonging to atom */
{
  gentop_qalpha *qa = (gentop_qalpha *) qat;
  int i;
  
  for(i=0; (i<qa->nr); i++)
    if (strcasecmp(atom,qa->qa[i].atom) == 0)
      return qa->qa[i].q;
      
  gmx_fatal(FARGS,"Can not find charge for atom %s",atom);
  
  return 0;
}

double get_qa_alpha(char *atom,gentop_qat qat)
/* Return the alpha belonging to atom */
{
  gentop_qalpha *qa = (gentop_qalpha *) qat;
  int i;
  
  for(i=0; (i<qa->nr); i++)
    if (strcasecmp(atom,qa->qa[i].atom) == 0)
      return qa->qa[i].alpha;
      
  gmx_fatal(FARGS,"Can not find polarizability for atom %s",atom);
  
  return 0;
}
