/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * GROtesk MACabre and Sinister
 */

#include <stdio.h>
#include <string.h>
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
  long fpos;

  fpos=ftell(fp);
  nblockwrite(fp,mtab->msize,mtab->mtab);
  return (ftell(fp)-fpos);
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

