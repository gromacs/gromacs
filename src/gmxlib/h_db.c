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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_h_db_c = "$Id$";

#include <string.h>
#include "string2.h"
#include "sysstuff.h"
#include "smalloc.h"
#include "futil.h"
#include "symtab.h"
#include "h_db.h"

/* There are 9 types of adding hydrogens, numbered from
 * 1 thru 9. Each of these has a specific number of
 * control atoms, that determine how the hydrogens are added.
 * Here these number are given. Because arrays start at 0 an
 * extra dummy for index 0 is added 
 */
/*                       1  2  3  4  5  6  7  8  9 */
int ncontrol[12] = { -1, 3, 3, 3, 3, 4, 3, 1, 3, 3 };

int compaddh(const void *a,const void *b)
{
  t_hackblock *ah,*bh;

  ah=(t_hackblock *)a;
  bh=(t_hackblock *)b;
  return strcasecmp(ah->name,bh->name);
}

void read_ab(char *line,char *fn,t_hack *hack)
{
  int  i,n,nh,tp,ncntl;
  char buf[80];
  
  if (sscanf(line,"%d%d%n",&nh,&tp,&n) != 2)
    fatal_error(0,"wrong format in input file %s on line\n%s\n",fn,line);
  line+=n;
  hack->nr=nh;
  hack->tp=tp;
  ncntl=ncontrol[tp];
  for(i=0; i<ncntl; i++) {
    if (sscanf(line,"%s%n",buf,&n) != 1)
      fatal_error(0,"Expected %d control atoms instead of %d while reading Hydrogen Database %s on line\n%s\n",ncntl,i-1,fn,line);
    hack->a[i]=strdup(buf);
    line+=n;
  }
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
  while (fgets(line,STRLEN,in)) {
    sscanf(line,"%s%n",buf,&n);
    if (debug) fprintf(debug,"%s",buf);
    srenew(aah,++nah);
    aah[nah-1].name=strdup(buf);
    sscanf(line+n,"%d",&nab);
    if (debug) fprintf(debug,"\t%d\n",nab);
    snew(aah[nah-1].hack,nab);
    aah[nah-1].nhack=nab;
    for(i=0; (i<nab); i++) {
      if (feof(in))
	fatal_error(0, "Expected %d lines of hydrogens, found only %d "
		    "while reading Hydrogen Database %s residue %s",
		    nab, i-1, aah[nah-1].name, hfn);
      fgets(buf, STRLEN, in);
      read_ab(buf,hfn,&(aah[nah-1].hack[i]));
      if (debug) print_ab(debug, &(aah[nah-1].hack[i]));
    }
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
  for(i=0; i < ncontrol[hack->tp]; i++)
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
