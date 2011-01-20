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

#include <stdio.h>
#include <string.h>
#include "sysstuff.h"
#include "string2.h"
#include "typedefs.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "txtdump.h"
#include "symtab.h"
#include "macros.h"

#define	BUFSIZE			1024
#define	TABLESIZE		5

static char *trim_string(const char *s,char *out, int maxlen)
     /*
      * Returns a pointer to a static area which contains a copy 
      * of s without leading or trailing spaces. Strings are
      * truncated to BUFSIZE positions.
      */      
{
  int len,i;
 
  if(strlen(s)>(size_t)(maxlen-1))
    gmx_fatal(FARGS,"Character buffer size too small\n");
  
  for (; (*s)&&((*s)==' '); s++);
  for (len=strlen(s); (len>0); len--) if (s[len-1]!=' ') break;
  if (len>=BUFSIZE) len=BUFSIZE-1;
  for (i=0; i<len; i++) out[i]=*(s++);
  out[i]=0;
  return out;
}

int lookup_symtab(t_symtab *symtab,char **name)
{
  int base,index;
  t_symbuf *symbuf;
  
  base=0;
  index=0;
  symbuf=symtab->symbuf;
  while (symbuf!=NULL) {
    index=name-symbuf->buf;
    if ( ( index >= 0 ) && ( index < symbuf->bufsize ) )
      return index+base;
    else {
      base+=symbuf->bufsize;
      symbuf=symbuf->next;
    }
  }
  gmx_fatal(FARGS,"symtab lookup \"%s\" not found",*name);
  return -1;
}

char **get_symtab_handle(t_symtab *symtab,int name)
{
  t_symbuf *symbuf;
  
  symbuf=symtab->symbuf;
  while (symbuf!=NULL) {
    if (name<symbuf->bufsize)
      return &(symbuf->buf[name]);
    else {
      name-=symbuf->bufsize;
      symbuf=symbuf->next;
    }
  }
  gmx_fatal(FARGS,"symtab get_symtab_handle %d not found",name);
  return NULL;
}

static t_symbuf *new_symbuf(void)
{
  t_symbuf *symbuf;

  snew(symbuf,1);
  symbuf->bufsize=TABLESIZE;
  snew(symbuf->buf,symbuf->bufsize);
  symbuf->next=NULL;

  return symbuf;
}

static char **enter_buf(t_symtab *symtab,char *name)
{
  int      i;
  t_symbuf *symbuf;
  gmx_bool     bCont;
  
  if (symtab->symbuf == NULL)
    symtab->symbuf=new_symbuf();

  symbuf=symtab->symbuf;
  do {
    for(i=0; (i < symbuf->bufsize); i++) {
      if (symbuf->buf[i]==NULL) {
	symtab->nr++;
	symbuf->buf[i]=strdup(name);
	return &(symbuf->buf[i]);
      } else if (strcmp(symbuf->buf[i],name)==0)
	return &(symbuf->buf[i]);
    }
    if (symbuf->next != NULL) {
      symbuf=symbuf->next;
      bCont = TRUE;
    }
    else
      bCont = FALSE;
  } while (bCont);

  symbuf->next=new_symbuf();
  symbuf=symbuf->next;

  symtab->nr++;
  symbuf->buf[0]=strdup(name);
  return &(symbuf->buf[0]);
}

char **put_symtab(t_symtab *symtab,const char *name)
{
  char buf[256];
  
  return enter_buf(symtab,trim_string(name,buf,255));
}

void open_symtab(t_symtab *symtab)
{
  symtab->nr=0;
  symtab->symbuf=NULL;
}

void close_symtab(t_symtab *symtab)
{
}

void done_symtab(t_symtab *symtab)
{
  int i;
  t_symbuf *symbuf,*freeptr;
  
  close_symtab(symtab);
  symbuf=symtab->symbuf;
  while (symbuf!=NULL) {
    for (i=0; (i < symbuf->bufsize) && (i < symtab->nr); i++)
      sfree(symbuf->buf[i]);
    symtab->nr-=i;
    sfree(symbuf->buf);
    freeptr=symbuf;
    symbuf=symbuf->next;
    sfree(freeptr);
  }
  symtab->symbuf=NULL;
  if (symtab->nr != 0)
    gmx_incons("Freeing symbol table (symtab) structure");
}

void free_symtab(t_symtab *symtab)
{
  t_symbuf *symbuf,*freeptr;
  
  close_symtab(symtab);
  symbuf=symtab->symbuf;
  while (symbuf!=NULL) {
    symtab->nr-=min(symbuf->bufsize,symtab->nr);
    freeptr=symbuf;
    symbuf=symbuf->next;
    sfree(freeptr);
  }
  symtab->symbuf=NULL;
  if (symtab->nr != 0)
    gmx_incons("Freeing symbol table (symtab) structure");
}

void pr_symtab(FILE *fp,int indent,const char *title,t_symtab *symtab)
{
  int i,j,nr;
  t_symbuf *symbuf;
  
  if (available(fp,symtab,indent,title)) 
    {
      indent=pr_title_n(fp,indent,title,symtab->nr);
      i=0;
      nr=symtab->nr;
      symbuf=symtab->symbuf;
      while (symbuf!=NULL)
        {
          for (j=0; (j < symbuf->bufsize) && (j < nr); j++)
            {
              pr_indent(fp,indent);
              (void) fprintf(fp,"%s[%d]=\"%s\"\n",title,i++,symbuf->buf[j]);
            }
          nr-=j;
          symbuf=symbuf->next;
        }
      if (nr != 0)
	gmx_incons("Printing symbol table (symtab) structure");
    }
}
