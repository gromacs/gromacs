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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_xlate_c = "$Id$";

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
  static char  *xlfile="xlateat.dat";
  
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

void rename_atoms(t_atoms *atoms,t_symtab *symtab)
{
  int nxlate,a,i;
  t_xlate_atom *xlatom;
  char c,*res,atombuf[32];
  bool bRenamed;

  xlatom = get_xlatoms(&nxlate);

  for(a=0; (a<atoms->nr); a++) {
    res  = *(atoms->resname[atoms->atom[a].resnr]);
    strcpy(atombuf,*(atoms->atomname[a]));
    if (isdigit(atombuf[0])) {
      c = atombuf[0];
      for (i=0; (i<strlen(atombuf)-1); i++)
	atombuf[i]=atombuf[i+1];
      atombuf[i]=c;
    }
    bRenamed=FALSE;
    for(i=0; (i<nxlate) && !bRenamed; i++) {
      if ((xlatom[i].res == NULL) || (strcasecmp(res,xlatom[i].res) == 0) ||
	  ((strcasecmp("protein",xlatom[i].res) == 0) && is_protein(res)))
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

