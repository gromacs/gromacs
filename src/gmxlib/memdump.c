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
          fprintf(log,"%.8X :",i);
          index=0;
        }
      (void) fprintf(log," %.2X",b);
      s[index++]=ascii(b);
    }
  print_chars(log,s,index);
  (void) fflush(log);
}
