/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
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
 * S  C  A  M  O  R  G
 */
static char *SRCID_xlate_c = "$Id$";

#include <ctype.h>
#include <string.h>
#include "typedefs.h"
#include "pdbio.h"
#include "strdb.h"
#include "string2.h"
#include "smalloc.h"
#include "index.h"
#include "xlate.h"

typedef struct {
  char *res;
  char *atom;
  char *replace;
} t_xlate_atom;

static t_xlate_atom *get_xlatoms(int *nxlatom)
{
  static char  *xlfile="xlateat.dat";
  
  t_xlate_atom *xl;
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

void xlate_atom(char *res,char **atom)
{
  static int          nxlate=0;
  static t_xlate_atom *xlatom;
  int  i;
  char c;
  
  if (nxlate == 0) {
    xlatom = get_xlatoms(&nxlate);
  }
  if (isdigit(*atom[0])) {
    c=*atom[0];
    for (i=0; (i<strlen(*atom)-1); i++)
      *atom[i]=*atom[i+1];
    srenew(*atom,strlen(*atom)+2);
    *atom[i]=c;
  }
  for(i=0; (i<nxlate); i++) {
    if ((xlatom[i].res == NULL) || (strcasecmp(res,xlatom[i].res) == 0) ||
	((strcasecmp("protein",xlatom[i].res) == 0) && is_protein(res)))
      if (strcasecmp(*atom,xlatom[i].atom) == 0) {
	sfree(*atom);
	*atom=strdup(xlatom[i].replace);
	return;
      }
  }
}

void rename_atoms(t_atoms *atoms)
{
  int i;
  
  for(i=0; (i<atoms->nr); i++) 
    xlate_atom(*atoms->resname[atoms->atom[i].resnr],atoms->atomname[i]);
}

