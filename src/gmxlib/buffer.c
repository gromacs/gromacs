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
 * Gromacs Runs On Most of All Computer Systems
 */
static char *SRCID_buffer_c = "$Id$";
#include "buffer.h"

extern int filler;

int mask(int i)
{
/*
  return filler;
*/
  return (i&0xff)|(((i+1)&0xff)<<8)|(((i+2)&0xff)<<16)|(((i+3)&0xff)<<24);
/*
  int result,j;

  i=-13-i;
  result=i<<12;
  i&=0x1f;
  for (j=0; j<i; j++) result=(result<<7)^i;
  return result;
*/
}

void clear_buff(int data[],int items)
{
  int i;

  for (i=0; i<items; i++) data[i]=0;
}

#define ERR_MASK 0xaa55aa55

void fill_buff(int data[],int items)
{
  int i;

  for (i=0; i<items; i++) data[i]=mask(i);
  /*for (i=0; i<items; i++) data[i]=ERR_MASK;*/
}

int check_buff(FILE *fp,char *title,int data[],int items,int verbose)
{
  int i,errs,check,comp;

  errs=0;
  if (!verbose)
    for (i=0; i<items; i++) { if (data[i]!=mask(i)) errs++; }
  else
    {
      for (i=0; i<items; i++) 
        {
          check=data[i];
          /* comp=ERR_MASK; */
          comp=mask(i);
          if (check!=comp)
            {
              fprintf(fp,"error: (%s) data: 0x%.8x, expected: 0x%.8x\n",
		      title,(unsigned int)check,(unsigned int)comp);
              errs++; 
            }
        }
    }
  return errs;
}
