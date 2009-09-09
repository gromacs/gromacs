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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <string.h>
#include "typedefs.h"
#include "strdb.h"
#include "string2.h"
#include "smalloc.h"
#include "symtab.h"
#include "index.h"

typedef struct {
  char *res;
  char *atom;
  char *replace;
} t_xlate_atom;

static t_xlate_atom *get_xlatoms(int *nxlatom)
{
  const char  *xlfile="xlateat.dat";
  
  t_xlate_atom *xl=NULL;
  char rbuf[32],abuf[32],repbuf[32];
  char **lines,*_ptr;
  int  nlines,i,n;
  
  nlines = get_lines(xlfile,&lines);
  if (nlines > 0) 
    snew(xl,nlines);
    
  n = 0;
  for(i=0; (i<nlines); i++) {
    if (sscanf(lines[i],"%s%s%s",rbuf,abuf,repbuf) != 3) 
      fprintf(stderr,"Invalid line '%s' in %s\n",lines[i],xlfile);
    else {
      /* Use wildcards... */
      if (strcmp(rbuf,"*") != 0)
	xl[n].res = strdup(rbuf);
      else
	xl[n].res = NULL;
      
      /* Replace underscores in the string by spaces */
      while ((_ptr = strchr(abuf,'_')) != 0)
	*_ptr = ' ';
      
      xl[n].atom = strdup(abuf);
      xl[n].replace = strdup(repbuf);
      n++;
    }
    sfree(lines[i]);
  }
  if (nlines > 0)
    sfree(lines);
  fprintf(stderr,"%d out of %d lines of %s converted succesfully\n",
	  n,nlines,xlfile);
  
  *nxlatom = n;
  
  return xl;
}

static void done_xlatom(int nxlate,t_xlate_atom **xlatom)
{
  int i;
  
  for(i=0; (i<nxlate); i++) {
    if ((*xlatom)[i].res)
      sfree((*xlatom)[i].res);
    if ((*xlatom)[i].atom)
      sfree((*xlatom)[i].atom);
    if ((*xlatom)[i].replace)
      sfree((*xlatom)[i].replace);
  }
  sfree(*xlatom);
  *xlatom = NULL;
}

void rename_atoms(t_atoms *atoms,t_symtab *symtab,t_aa_names *aan)
{
  int nxlate,a,i;
  t_xlate_atom *xlatom;
  char c,*res,atombuf[32];
  bool bRenamed;

  xlatom = get_xlatoms(&nxlate);

  for(a=0; (a<atoms->nr); a++) {
    res = *(atoms->resinfo[atoms->atom[a].resind].name);
    strcpy(atombuf,*(atoms->atomname[a]));
    if (isdigit(atombuf[0])) {
      c = atombuf[0];
      for (i=0; ((size_t)i<strlen(atombuf)-1); i++)
	atombuf[i]=atombuf[i+1];
      atombuf[i]=c;
    }
    bRenamed=FALSE;
    for(i=0; (i<nxlate) && !bRenamed; i++) {
      if ((xlatom[i].res == NULL) || (strcasecmp(res,xlatom[i].res) == 0) ||
	  ((strcasecmp("protein",xlatom[i].res) == 0) && is_protein(aan,res)))
	if (strcasecmp(atombuf,xlatom[i].atom) == 0) {
	  /* don't free the old atomname, since it might be in the symtab */
	  strcpy(atombuf,xlatom[i].replace);
	  bRenamed=TRUE;
	}
    }
    if (strcmp(atombuf,*atoms->atomname[a]) != 0)
      atoms->atomname[a] = put_symtab(symtab,atombuf);
  }

  done_xlatom(nxlate,&xlatom);
}

