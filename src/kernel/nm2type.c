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
#include "toppush.h"
#include "pdb2top.h"
#include "topexcl.h"
#include "x2top.h"

t_nm2type *rd_nm2type(char *ff,int *nnm)
{
  FILE      *fp;
  bool      bCont;
  char      libfilename[128];
  char      format[128],f1[128];
  char      buf[1024],elem[16],type[16],nbbuf[16],**newbuf;
  int       i,nb,nnnm,line=1;
  t_nm2type *nm2t=NULL;
  
  sprintf(libfilename,"%s.n2t",ff);
  fp = libopen(libfilename);
  
  nnnm = 0;
  do {
    /* Read a line from the file */
    bCont = (fgets2(buf,1023,fp) != NULL); 
    
    if (bCont) {
      /* Remove comment */
      strip_comment(buf);
      if (sscanf(buf,"%s%s%d",elem,type,&nb) == 3) {
	/* If we can read the first three, there probably is more */
	if (nb > 0) {
	  snew(newbuf,nb);
	  strcpy(format,"%*s%*s%*d");
	  for(i=0; (i<nb); i++) {
	    /* Complicated format statement */
	    strcpy(f1,format);
	    strcat(f1,"%s");
	    if (sscanf(buf,f1,nbbuf) != 1)
	      gmx_fatal(FARGS,"Error on line %d of %s",line,libfilename);
	    newbuf[i] = strdup(nbbuf);
	    strcat(format,"%*s");
	  }
	}
	else
	  newbuf = NULL;
	srenew(nm2t,nnnm+1);
	nm2t[nnnm].elem   = strdup(elem);
	nm2t[nnnm].type   = strdup(type);
	nm2t[nnnm].nbonds = nb;
	nm2t[nnnm].bond   = newbuf;
	nnnm++;
      }
      line++;
    }
  } while(bCont);
  fclose(fp);
  
  *nnm = nnnm;
  
  return nm2t;
}

void dump_nm2type(FILE *fp,int nnm,t_nm2type nm2t[])
{
  int i,j;
  
  fprintf(fp,"; nm2type database\n");
  for(i=0; (i<nnm); i++) {
    fprintf(fp,"%-8s%-8s%-4d",nm2t[i].elem,nm2t[i].type,nm2t[i].nbonds);
    for(j=0; (j<nm2t[i].nbonds); j++)
      fprintf(fp,"%-5s",nm2t[i].bond[j]);
    fprintf(fp,"\n");
  }
}

char *nm2type(int nnm,t_nm2type nm2t[],char *nm,int nbonds)
{
  int i;

  /* First check for names */  
  for(i=0; (i<nnm); i++) {
    if ((strcasecmp(nm2t[i].elem,nm) == 0) &&
	(nm2t[i].nbonds == nbonds))
      return nm2t[i].type;
  }
  /* Then for element */
  for(i=0; (i<nnm); i++) {
    if ((nm2t[i].elem[0] == nm[0]) &&
	(nm2t[i].nbonds == nbonds))
      return nm2t[i].type;
  }
	      
  return NULL;
}
     
