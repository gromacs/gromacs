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
 * Gnomes, ROck Monsters And Chili Sauce
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
  t_addh *ah,*bh;

  ah=(t_addh *)a;
  bh=(t_addh *)b;
  return strcasecmp(ah->resname,bh->resname);
}

t_addh *search_hb(char *key,int nh,t_addh ah[])
{
  t_addh ahkey;

  ahkey.resname=key;
  return (t_addh *)bsearch((void *)&ahkey,ah,nh,
			   (size_t)sizeof(ahkey),compaddh);
}

t_addh *search_ter_hb(char *key,int nh,t_addh ah[])
{
  t_addh ahkey;

  ahkey.resname=key;
  return (t_addh *)bsearch((void *)&ahkey,ah,nh,
			   (size_t)sizeof(ahkey),compaddh);
}

void read_ab(FILE *in,char *fn,t_add_block *ab)
{
  int  i,nh,tp,ncntl;
  char buf[80];

  if (fscanf(in,"%d%d",&nh,&tp) != 2)
    fatal_error(0,"error in input file %s\n",fn);
  ab->nh=nh;
  ab->tp=tp;
  ncntl=ncontrol[tp];
  if (debug) 
    fprintf(debug,"%s, %d: will read %d elements for type %d\n",
	    __FILE__,__LINE__,ncntl,tp);
  for(i=0; (i<ncntl); i++) {
    if (fscanf(in,"%s",buf) != 1) 
      fatal_error(0,"%d instead of %d (ncntl) elements found in input file\n",i,ncntl);
    ab->na[i]=strdup(buf);
  }
}

int read_h_db(char *fn,t_addh **ah)
{	
  FILE   *in;
  char   hfn[256];
  char   buf[80];
  int    i,nab,nah=0;
  t_addh *aah=NULL;

  sprintf(hfn,"%s.hdb",fn);
  in=libopen(hfn);
  while (fscanf(in,"%s",buf) == 1) {
    srenew(aah,++nah);
    aah[nah-1].resname=strdup(buf);
    fscanf(in,"%d",&nab);
    snew(aah[nah-1].ab,nab);
    aah[nah-1].n_add=nab;
    for(i=0; (i<nab); i++) 
      read_ab(in,hfn,&(aah[nah-1].ab[i]));
  }
  fclose(in);
  
  /* Sort the list (necessary to be able to use bsearch */
  qsort(aah,nah,(size_t)sizeof(**ah),compaddh);

  *ah=aah;
  return nah;
}

void print_ab(FILE *out,t_add_block *ab)
{
  int i,ncntl;

  fprintf(out,"\t%d\t%d",ab->nh,ab->tp);
  ncntl=ncontrol[ab->tp];
  for(i=0; (i<ncntl); i++)
    fprintf(out,"\t%s",ab->na[i]);
  fprintf(out,"\n");
}

void print_h_db(FILE *out,int nh,t_addh ah[])
{
  int i,j;

  for(i=0; (i<nh); i++) {
    fprintf(out,"%s\t%d\n",ah[i].resname,ah[i].n_add);
    for(j=0; (j<ah[i].n_add); j++)
      print_ab(out,&(ah[i].ab[j]));
  }
}

t_addh *search_h_db(int nh,t_addh ah[],char *key)
{
  t_addh ahkey,*result;

  if (nh <= 0)
    return NULL;
    
  ahkey.resname=key;

  result=(t_addh *)bsearch(&ahkey,ah,nh,(size_t)sizeof(ah[0]),compaddh);
  
  return result;
}
