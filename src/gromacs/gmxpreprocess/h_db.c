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
#include "gmxfio.h"
#include "fflibutil.h"
#include "gmx_fatal.h"

/* There are 11 types of adding hydrogens, numbered from
 * 1 thru 11. Each of these has a specific number of
 * control atoms, that determine how the hydrogens are added.
 * Here these number are given. Because arrays start at 0 an
 * extra dummy for index 0 is added 
 */
const int ncontrol[] = { -1, 3, 3, 3, 3, 4, 3, 1, 3, 3, 1, 1 };
#define maxcontrol asize(ncontrol)

int compaddh(const void *a,const void *b)
{
  t_hackblock *ah,*bh;

  ah=(t_hackblock *)a;
  bh=(t_hackblock *)b;
  return gmx_strcasecmp(ah->name,bh->name);
}

void print_ab(FILE *out,t_hack *hack,char *nname)
{
  int i;

  fprintf(out,"%d\t%d\t%s",hack->nr,hack->tp,nname);
  for(i=0; (i < hack->nctl); i++)
    fprintf(out,"\t%s",hack->a[i]);
  fprintf(out,"\n");
}


void read_ab(char *line,const char *fn,t_hack *hack)
{
  int  i,nh,tp,ns;
  char a[4][12];
  char hn[32];
  
  ns = sscanf(line,"%d%d%s%s%s%s%s",&nh,&tp,hn,a[0],a[1],a[2],a[3]);
  if (ns < 4)
    gmx_fatal(FARGS,"wrong format in input file %s on line\n%s\n",fn,line);
  
  hack->nr=nh;
  hack->tp=tp;
  if ((tp < 1) || (tp >= maxcontrol)) 
    gmx_fatal(FARGS,"Error in hdb file %s:\nH-type should be in 1-%d. Offending line:\n%s",fn,maxcontrol-1,line);
  
  hack->nctl = ns - 3;
  if ((hack->nctl != ncontrol[hack->tp]) && (ncontrol[hack->tp] != -1))
    gmx_fatal(FARGS,"Error in hdb file %s:\nWrong number of control atoms (%d iso %d) on line:\n%s\n",fn,hack->nctl,ncontrol[hack->tp],line);
  for(i=0; (i<hack->nctl); i++) {
    hack->a[i]=strdup(a[i]);
  }
  for(   ; i<4; i++) {
    hack->a[i]=NULL;
  }
  hack->oname=NULL;
  hack->nname=strdup(hn);
  hack->atom=NULL;
  hack->cgnr=NOTSET;
  hack->bXSet=FALSE;
  for(i=0; i<DIM; i++)
    hack->newx[i]=NOTSET;
}

static void dump_h_db(const char *fn,int nah,t_hackblock *ah)
{
  FILE *fp;
  char buf[STRLEN],nname[STRLEN];
  int  i,j,k;
  
  sprintf(buf,"%s_new.hdb",fn);
  fp = gmx_fio_fopen(buf,"w");
  for(i=0; (i<nah); i++) {
    fprintf(fp,"%-8s%-8d\n",ah[i].name,ah[i].nhack);
    for(k=0; (k<ah[i].nhack); k++) {
      strcpy(nname,ah[i].hack[k].a[0]);
      nname[0] = 'H';
      print_ab(fp,&ah[i].hack[k],nname);
    }
  }
  gmx_fio_fclose(fp);
}

static void read_h_db_file(const char *hfn,int *nahptr,t_hackblock **ah)
{	
  FILE   *in;
  char   filebase[STRLEN],line[STRLEN], buf[STRLEN];
  int    i, n, nab, nah;
  t_hackblock *aah;

  if (debug) fprintf(debug,"Hydrogen Database (%s):\n",hfn);

  fflib_filename_base(hfn,filebase,STRLEN);
  /* Currently filebase is read and set, but not used.
   * hdb entries from any hdb file and be applied to rtp entries
   * in any rtp file.
   */

  in = fflib_open(hfn);

  nah = *nahptr;
  aah = *ah;
  while (fgets2(line,STRLEN-1,in)) {
    if (sscanf(line,"%s%n",buf,&n) != 1) {
      fprintf(stderr,"Error in hdb file: nah = %d\nline = '%s'\n",
	      nah,line);
      break;
    }
    if (debug) fprintf(debug,"%s",buf);
    srenew(aah,nah+1);
    clear_t_hackblock(&aah[nah]);
    aah[nah].name     = strdup(buf);
    aah[nah].filebase = strdup(filebase);
    
    if (sscanf(line+n,"%d",&nab) == 1) {
      if (debug) fprintf(debug,"  %d\n",nab);
      snew(aah[nah].hack,nab);
      aah[nah].nhack = nab;
      for(i=0; (i<nab); i++) {
	if (feof(in))
	  gmx_fatal(FARGS, "Expected %d lines of hydrogens, found only %d "
		      "while reading Hydrogen Database %s residue %s",
		      nab, i-1, aah[nah].name, hfn);
	if(NULL==fgets(buf, STRLEN, in))
        {
	  gmx_fatal(FARGS,"Error reading from file %s",hfn);
	}
	read_ab(buf,hfn,&(aah[nah].hack[i]));
      }
    }
    nah++;
  }
  ffclose(in);
  
  /* Sort the list (necessary to be able to use bsearch */
  qsort(aah,nah,(size_t)sizeof(**ah),compaddh);

  if (debug)
    dump_h_db(hfn,nah,aah);
  
  *nahptr = nah;
  *ah     = aah;
}

int read_h_db(const char *ffdir,t_hackblock **ah)
{
  int  nhdbf,f;
  char **hdbf;
  int  nah;
  FILE *fp;

  /* Read the hydrogen database file(s).
   * Do not generate an error when no files are found.
   */
  nhdbf = fflib_search_file_end(ffdir,".hdb",FALSE,&hdbf);
  nah = 0;
  *ah = NULL;
  for(f=0; f<nhdbf; f++) {
    read_h_db_file(hdbf[f],&nah,ah);
    sfree(hdbf[f]);
  }
  sfree(hdbf);

  return nah;
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
