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
