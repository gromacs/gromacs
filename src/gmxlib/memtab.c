/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Glycine aRginine prOline Methionine Alanine Cystine Serine
 */
static char *SRCID_memtab_c = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "sysstuff.h"
#include "fatal.h"
#include "smalloc.h"
#include "memtab.h"
#include "binio.h"
        
#define	ALIGN_BASE	sizeof(double)
#define	ALIGN(size)	(((ALIGN_BASE-1+(size))/ALIGN_BASE)*ALIGN_BASE)

static char *extend_memtab(t_memtab *mtab,int size,void *ref,void **handle)
{
  int i,offset,aligned_size;
  char *newpart;

  i=mtab->nref++;
  offset=mtab->msize;
  aligned_size=ALIGN(size);
  mtab->refs[i].ref=ref;
  mtab->refs[i].handle=(void *)offset;
  *handle=mtab->refs[i].handle;
  mtab->msize+=aligned_size;
  srenew(mtab->mtab,mtab->msize);
  newpart=&mtab->mtab[offset];
  memset(newpart,0,aligned_size);
  return newpart;
}

void init_memtab(t_memtab *mtab)
{
  mtab->nref=0;
  mtab->refs=NULL;
  mtab->msize=0;
  mtab->mtab=NULL;
}

void dispose_memtab(t_memtab *mtab)
{
  sfree(mtab->refs);
  sfree(mtab->mtab);
  init_memtab(mtab);
}

void *put_memtab(t_memtab *mtab,int size,void *ref)
{
  int i;
  void *handle;
  
  if ((ref==NULL)||(size<0))
    return (void *)NOENTRY;
  else
    {
      for (i=0; i<mtab->nref; i++)
        if (mtab->refs[i].ref==ref) return mtab->refs[i].handle;
      (void) memcpy(extend_memtab(mtab,size,ref,&handle),ref,size);
      return handle;
    }
}

void *get_memtab(t_memtab *mtab,void *handle)
{
  if (handle==(void *)NOENTRY)
    return NULL;
  else 
    return &(mtab->mtab[(int)handle]);
}

long wr_memtab(FILE *fp,t_memtab *mtab)
{
  off_t fpos;

#ifdef HAVE_FSEEKO  
  fpos=ftello(fp);
#else
  fpos=ftell(fp);
#else  
  nblockwrite(fp,mtab->msize,mtab->mtab);
#ifdef HAVE_FSEEKO  
  return (ftello(fp)-fpos);
#else
  return (ftell(fp)-fpos);
#endif  
}

void *rd_memtab(FILE *fp,int size,t_memtab *mtab)
{
  char *ref;
  void *handle;
  
  if (size<0)
    return (void *)NOENTRY;
  else
    {
      ref=extend_memtab(mtab,size,NULL,&handle);
      nblockread(fp,size,ref);
      return handle;
    }
}

