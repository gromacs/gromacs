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

#include <string.h>
#include "string2.h"
#include "sysstuff.h"
#include "smalloc.h"
#include "futil.h"
#include "symtab.h"
#include "h_db.h"

/* There are 11 types of adding hydrogens, numbered from
 * 1 thru 11. Each of these has a specific number of
 * control atoms, that determine how the hydrogens are added.
 * Here these number are given. Because arrays start at 0 an
 * extra dummy for index 0 is added 
 */
/* const int ncontrol[12] = { -1, 3, 3, 3, 3, 4, 3, 1, 3, 3, 1, 1 }; */

int compaddh(const void *a,const void *b)
{
  t_hackblock *ah,*bh;

  ah=(t_hackblock *)a;
  bh=(t_hackblock *)b;
  return strcasecmp(ah->name,bh->name);
}

void read_ab(char *line,char *fn,t_hack *hack)
{
  int  i,nh,tp,ns;
  char a[4][12];
  
  ns = sscanf(line,"%d%d%s%s%s%s",&nh,&tp,a[0],a[1],a[2],a[3]);
  if (ns < 3)
    gmx_fatal(FARGS,"wrong format in input file %s on line\n%s\n",fn,line);
  
  hack->nr=nh;
  hack->tp=tp;
  hack->nctl = ns - 2;
  for(i=0; (i<hack->nctl); i++) 
    hack->a[i]=strdup(a[i]);
  for(   ; i<4; i++)
    hack->a[i]=NULL;
  hack->oname=NULL;
  hack->nname=NULL;
  hack->atom=NULL;
  hack->cgnr=NOTSET;
  for(i=0; i<DIM; i++)
    hack->newx[i]=NOTSET;
}

int read_h_db(char *fn,t_hackblock **ah)
{	
  FILE   *in;
  char   hfn[STRLEN], line[STRLEN], buf[STRLEN];
  int    i, n, nab, nah;
  t_hackblock *aah;

  sprintf(hfn,"%s.hdb",fn);
  in=libopen(hfn);
  if (debug) fprintf(debug,"Hydrogen Database (%s):\n",hfn);
  nah=0;
  aah=NULL;
  while (fgets2(line,STRLEN-1,in)) {
    if (sscanf(line,"%s%n",buf,&n) != 1) {
      fprintf(stderr,"Error in hdb file: nah = %d\nline = '%s'\n",
	      nah,line);
      break;
    }
    if (debug) fprintf(debug,"%s",buf);
    srenew(aah,nah+1);
    clear_t_hackblock(&aah[nah]);
    aah[nah].name=strdup(buf);
    
    if (sscanf(line+n,"%d",&nab) == 1) {
      if (debug) fprintf(debug,"  %d\n",nab);
      snew(aah[nah].hack,nab);
      aah[nah].nhack = nab;
      for(i=0; (i<nab); i++) {
	if (feof(in))
	  gmx_fatal(FARGS, "Expected %d lines of hydrogens, found only %d "
		      "while reading Hydrogen Database %s residue %s",
		      nab, i-1, aah[nah].name, hfn);
	fgets(buf, STRLEN, in);
	read_ab(buf,hfn,&(aah[nah].hack[i]));
	if (debug) print_ab(debug, &(aah[nah].hack[i]));
      }
    }
    nah++;
  }
  fclose(in);
  
  /* Sort the list (necessary to be able to use bsearch */
  qsort(aah,nah,(size_t)sizeof(**ah),compaddh);

  *ah=aah;
  return nah;
}

void print_ab(FILE *out,t_hack *hack)
{
  int i;

  fprintf(out,"%d\t%d",hack->nr,hack->tp);
  for(i=0; (i < hack->nctl); i++)
    fprintf(out,"\t%s",hack->a[i]);
  fprintf(out,"\n");
}

void print_h_db(FILE *out,int nh,t_hackblock ah[])
{
  int i,j;

  for(i=0; (i<nh); i++) {
    fprintf(out,"%s\t%d\n",ah[i].name,ah[i].nhack);
    for(j=0; (j<ah[i].nhack); j++) {
      fprintf(out,"\t");
      print_ab(out,&(ah[i].hack[j]));
    }
  }
}

t_hackblock *search_h_db(int nh,t_hackblock ah[],char *key)
{
  t_hackblock ahkey,*result;

  if (nh <= 0)
    return NULL;
  
  ahkey.name=key;

  result=(t_hackblock *)bsearch(&ahkey,ah,nh,(size_t)sizeof(ah[0]),compaddh);
  
  return result;
}
