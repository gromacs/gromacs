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
static char *SRCID_smalloc_c = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fatal.h"
#include "smalloc.h"
#include "main.h"

#ifdef DEBUG
#define NN "NULL"
static void log_action(int bMal,char *what,char *file,int line,
                       int nelem,int size,void *ptr)
{
  static int btot=0;
  int        bytes;
  
  bytes=size*nelem;
  if (!bMal)
    bytes=-bytes;
  btot+=bytes;
    
  bytes/=1024;
#ifdef _amb_
  if ((stdlog != NULL) && ((bytes != 0) || ((unsigned long) ptr == 0x5490) )) {
    long mm=0,mx=0;
    
    mm=memavail();
    mx=maxavail();
    
    fprintf(stdlog,"%30s:%4d b (mm:%4d, mx:%4d) [%s, line %d, nelem %d, size %d, ptr: %x]\n",
	    what ? what : NN,bytes,mm/1024,mx/1024,
	    file ? file : NN,line,nelem,size,ptr);
  }
#else
  if ((stdlog != NULL) && (bytes != 0))
    fprintf(stdlog,"%30s:%6d kb (%7d kb) [%s, line %d, nelem %d, size %d]\n",
	    what ? what : NN,bytes,btot/1024,
	    file ? file : NN,line,nelem,size);
#endif
}
#undef NN
#endif

void *save_malloc(char *name,char *file,int line,int size)
{
  void *p;
  
  p=NULL;
  if (size==0)
    p=NULL;
  else
    {
      if ((p=malloc(size))==NULL) 
        fatal_error(errno,"malloc for %s (%d bytes, file %s, line %d)",
                    name,size,file,line);
      (void) memset(p,0,size);
    }
#ifdef DEBUG
  log_action(1,name,file,line,1,size,p);
#endif
  return p;
}

void *save_calloc(char *name,char *file,int line,
                  unsigned nelem,unsigned elsize)
{
  void *p;
  
  p=NULL;
  if ((nelem==0)||(elsize==0))
    p=NULL;
  else
    {
      if ((p=calloc((size_t)nelem,(size_t)elsize))==NULL) 
        fatal_error(errno,"calloc for %s (nelem=%d, elsize=%d, file %s"
                    ", line %d)",name,nelem,elsize,file,line);
    }
#ifdef DEBUG
  log_action(1,name,file,line,nelem,elsize,p);
#endif
  return p;
}

void *save_realloc(char *name,char *file,int line,void *ptr,unsigned size)
{
  void *p;
  
  p=NULL;
  if (size==0)
    p=NULL;
  else
    {
      if (ptr==NULL) 
	p=malloc((size_t)size); 
      else 
	p=realloc(ptr,(size_t)size);
      if (p==NULL) 
        fatal_error(errno,
                    "realloc for %s (%d bytes, file %s, line %d, %s=0x%8x)",
                    name,size,file,line,name,ptr);
    }
#ifdef DEBUG
  log_action(1,name,file,line,1,size,p);
#endif
  return p;
}

void save_free(char *name,char *file,int line,void *ptr)
{
#ifdef DEBUG
  log_action(0,name,file,line,0,0,ptr);
#endif
  if (ptr!=NULL)
    {
#ifdef _sun_
     if (free(ptr)==-1)
       fatal_error(errno,"free for %s (file %s, line %d, %s=0x%8x)",
                   name,file,line,name,ptr);
#else
     free(ptr);
#endif
   }
}

unsigned maxavail(void)
{
  char *ptr;
  unsigned low,high,size;
  
  low=0;
  high=256e6;
  while ((high-low)>4)
    {
      size=(high+low)/2;
      if ((ptr=malloc((size_t)size))==NULL)
        high=size;
      else
        {
          free(ptr);
          low=size;
        }
    }
  return low;
}

unsigned memavail(void)
{
  char *ptr;
  unsigned size;
  
  size=maxavail();
  if (size!=0)
    {
      if ((ptr=malloc((size_t)size))!=NULL)
        {
          size+=memavail();
          free(ptr);
        }
    }
  return size;
}
