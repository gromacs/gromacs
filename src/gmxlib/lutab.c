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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
static char *SRCID_lutab_c = "$Id$";

/*
 *	Invsqrt lookup tables @(#) gentab.c 1.9 11/5/92
 *
 *	Copyright (c) 1992, A. Sijbers, Groningen University
 */
#include <math.h>
#include "lutab.h"
#include "callf77.h"

/* Global variable that is exported */
t_lutab lookup_table;

static void init_table(t_lutab *lookup_table)
{
  word index,addr;
  int new_exp,i,exp;
  t_convert bit_in,bit_out;
  
  for (i=0; i<EXP_SEED_SIZE; i++)
    lookup_table->exp_seed[i]=NOT_INITIALISED;
  for (i=0; i<FRACT_SEED_SIZE; i++) 
    lookup_table->fract_seed[i]=NOT_INITIALISED;
  for (i=0; i<EXP_SEED_SIZE; i++)
    {
      if ((exp=(i-127))<0) new_exp=127-((exp+1)/2); else new_exp=127-(exp/2)-1;
      lookup_table->exp_seed[i]=((new_exp)<<EXP_SHIFT)&EXP_MASK;
    }
  index=FRACT_FIRST; 
  for (i=0; i<FRACT_SEED_SIZE; i++)
    {
      bit_in.bval=(index<<FRACT_SHIFT);
      addr=FRACT_ADDR(bit_in.bval);
      bit_out.fval=(1.0/sqrt(bit_in.fval));
      lookup_table->fract_seed[addr]=(bit_out.bval&FRACT_MASK);
      if (lookup_table->fract_seed[addr]==0) 
        lookup_table->fract_seed[addr]=(MAX_FRACT&FRACT_MASK);
      index++;
    }
}

void init_lookup_table(FILE *log)
{
#ifdef CINVSQRT
  fprintf(log,"Generating lookup table for invsqrt calculation in C\n");
  init_table(&lookup_table);
#endif
#ifdef USEF77
#ifdef FINVSQRT
  fprintf(log,"Generating lookup table for invsqrt calculation in Fortran\n");
  fillbuf();
#endif
#endif
}

