/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
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
 * GROup of MAchos and Cynical Suckers
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
              (void) fprintf(fp,"error: (%s) data: 0x%.8x, "
                             "expected: 0x%.8x\n",title,check,comp);
              errs++; 
            }
        }
    }
  return errs;
}
