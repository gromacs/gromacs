/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Gromacs Runs On Most of All Computer Systems
 */
static char *SRCID_symtab_c = "$Id$";

#include <stdio.h>
#include <string.h>
#include "sysstuff.h"
#include "string2.h"
#include "binio.h"
#include "assert.h"
#include "typedefs.h"
#include "fatal.h"
#include "smalloc.h"
#include "txtdump.h"
#include "symtab.h"

#define	BUFSIZE			1024
#define	TABLESIZE		5

static char *trim_string(char *s)
     /*
      * Returns a pointer to a static area which contains a copy 
      * of s without leading or trailing spaces. Strings are
      * truncated to BUFSIZE positions.
      */      
{
  static char buf[BUFSIZE];
  int len,i;
  
  for (; (*s)&&((*s)==' '); s++);
  for (len=strlen(s); (len>0); len--) if (s[len-1]!=' ') break;
  if (len>=BUFSIZE) len=BUFSIZE-1;
  for (i=0; i<len; i++) buf[i]=*(s++);
  buf[i]=0;
  return buf;
}

int lookup_symtab(t_symtab *symtab,char **name)
{
  int base,index;
  t_symbuf *symbuf;
  
  base=0;
  index=0;
  symbuf=symtab->symbuf;
  while (symbuf!=NULL)
    {
      index=name-symbuf->buf;
      if ((index>=0)&&(index<symbuf->bufsize))
        return index+base;
      else
        {
          base+=symbuf->bufsize;
          symbuf=symbuf->next;
        }
    }
  fatal_error(0,"symtab lookup \"%s\" not found",*name);
  return -1;
}

char **get_symtab_handle(t_symtab *symtab,int name)
{
  t_symbuf *symbuf;
  
  symbuf=symtab->symbuf;
  while (symbuf!=NULL)
    {
      if (name<symbuf->bufsize)
        return &(symbuf->buf[name]);
      else
        {
          name-=symbuf->bufsize;
          symbuf=symbuf->next;
        }
    }
  fatal_error(0,"symtab get_symtab_handle %d not found",name);
  return NULL;
}

static t_symbuf *new_symbuf(void)
{
  t_symbuf *dummy;

  snew(dummy,1);
  dummy->bufsize=TABLESIZE;
  snew(dummy->buf,dummy->bufsize); /* buf[i]==NULL ! */
  dummy->next=NULL;

  return dummy;
}

static char **enter_buf(t_symtab *symtab,char *name)
{
  int      i;
  t_symbuf *dummy;
  
  if (symtab->symbuf == NULL)
    symtab->symbuf=new_symbuf();

  dummy=symtab->symbuf;
  do {
    for(i=0; (i<dummy->bufsize); i++) {
      if (dummy->buf[i]==NULL) {
	symtab->nr++;
	dummy->buf[i]=strdup(name);
	return &(dummy->buf[i]);
      }
      else if (strcmp(dummy->buf[i],name)==0)
	return &(dummy->buf[i]);
    }
    if (dummy->next != NULL)
      dummy=dummy->next;
    else
      break;
  } while (1);

  dummy->next=new_symbuf();
  dummy=dummy->next;

  symtab->nr++;
  dummy->buf[0]=strdup(name);
  return &(dummy->buf[0]);
}

char **put_symtab(t_symtab *symtab,char *name)
{
  return enter_buf(symtab,trim_string(name));
}

void open_symtab(t_symtab *symtab)
{
  symtab->nr=0;
  symtab->symbuf=NULL;
}

void close_symtab(t_symtab *symtab)
{
}

void rm_symtab(t_symtab *symtab)
{
  int i;
  t_symbuf *symbuf,*freeptr;
  
  close_symtab(symtab);
  symbuf=symtab->symbuf;
  while (symbuf!=NULL)
    {
      for (i=0; (i<symbuf->bufsize)&&(i<symtab->nr); i++)
        sfree(symbuf->buf[i]);
      symtab->nr-=i;
      sfree(symbuf->buf);
      freeptr=symbuf;
      symbuf=symbuf->next;
      sfree(freeptr);
    }
  symtab->symbuf=NULL;
  assert(symtab->nr==0);
}

long wr_symtab(FILE *fp,t_symtab *symtab)
{
  int i,nr,len;
  long fpos;
  t_symbuf *symbuf;

  fpos=ftell(fp);
  blockwrite(fp,symtab->nr);
  nr=symtab->nr;
  symbuf=symtab->symbuf;
  while (symbuf!=NULL)
    {
      for (i=0; (i<symbuf->bufsize)&&(i<nr); i++)
        {
          len=strlen(symbuf->buf[i])+1;
          blockwrite(fp,len);
          nblockwrite(fp,len,symbuf->buf[i]);
        }
      nr-=i;
      symbuf=symbuf->next;
    }
  assert(nr==0);
  return (ftell(fp)-fpos);
}

long rd_symtab(FILE *fp,t_symtab *symtab)
{
  int i,nr,len;
  long fpos;

  fpos=ftell(fp);
  blockread(fp,symtab->nr);
  nr=symtab->nr;
  snew(symtab->symbuf,1);
  symtab->symbuf->bufsize=nr;
  snew(symtab->symbuf->buf,nr);
  for (i=0; i<nr; i++)
    {
      blockread(fp,len);
      snew(symtab->symbuf->buf[i],len);
      nblockread(fp,len,symtab->symbuf->buf[i]);
    }
  return (ftell(fp)-fpos);
}

void pr_symtab(FILE *fp,int indent,char *title,t_symtab *symtab)
{
  int i,j,nr;
  t_symbuf *symbuf;
  
  if (available(fp,symtab,title))
    {
      indent=pr_title_n(fp,indent,title,symtab->nr);
      i=0;
      nr=symtab->nr;
      symbuf=symtab->symbuf;
      while (symbuf!=NULL)
        {
          for (j=0; (j<symbuf->bufsize)&&(j<nr); j++)
            {
              pr_indent(fp,indent);
              (void) fprintf(fp,"%s[%d]=\"%s\"\n",title,i++,symbuf->buf[j]);
            }
          nr-=j;
          symbuf=symbuf->next;
        }
      assert(nr==0);
    }
}
