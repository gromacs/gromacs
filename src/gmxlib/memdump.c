/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Gyas ROwers Mature At Cryogenic Speed
 */
static char *SRCID_memdump_c = "$Id$";
#include <stdio.h>
#include "memdump.h"

#define NUL                  '\0'
#define LINE_WIDTH           16
#define byte                 unsigned char

static void print_chars(FILE *fp,char s[],int len)
{
  int i;
  
  if (len)
    {
      for (i=len; i<LINE_WIDTH; i++) (void) fprintf(fp,"   ");
      s[len]=NUL;
      (void) fprintf(fp," | %s\n",s);
    }
}

static char ascii(byte b)
{     
  if ((b<' ')||(b>'~')) return ('.'); else return (char) b;
}

void mem_dump(FILE *log,char *title,void *mem,int len)
{
  int i,index;
  char s[LINE_WIDTH+1];
  byte *p,b;
  
  p=mem;
  index=0;
  s[0]=NUL;
  if (len) (void) fprintf(log,"memdump of %s:\n",title);
  for (i=0; i<len; i++)
    {
      b=p[i];
      if ((index==0)||(index==LINE_WIDTH))
        {
          print_chars(log,s,index); 
          fprintf(log,"%.8X :",(unsigned int)i);
          index=0;
        }
      (void) fprintf(log," %.2X",(unsigned int)b);
      s[index++]=ascii(b);
    }
  print_chars(log,s,index);
  (void) fflush(log);
}
