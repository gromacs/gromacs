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
 * GROningen MAchine for Chemical Simulation
 */
#include "sysstuff.h"
#include "typedefs.h"
#include "splittop.h"
#include "smalloc.h"
#include "fatal.h"
#include "main.h"

static void split_ilist(FILE *log,t_ilist *il,t_commrec *cr)
{
  t_iatom *ia;
  int     i,start,end,nr;
  
  if (cr->pid == 0)
    start=0;
  else
    start=il->multinr[cr->pid-1];
  end=il->multinr[cr->pid];
  
  nr=end-start;
  snew(ia,nr);

  for(i=0; (i<nr); i++)
    ia[i]=il->iatoms[start+i];

  sfree(il->iatoms);
  il->iatoms=ia;
  
  for(i=0; (i<MAXPROC); i++)
    il->multinr[i]=nr;
  il->nr=nr;
}

static void split_idef(FILE *log,t_idef *idef,t_commrec *cr)
{
  int i;
  
  for(i=0; (i<F_NRE); i++)
    split_ilist(log,&idef->il[i],cr);
}
	
void mdsplit_top(FILE *log,t_topology *top,t_commrec *cr)
{
  if (cr->nprocs < 2)
    return;
    
  split_idef(log,&top->idef,cr);
#ifdef DEBUG
  pr_idef(log,0,"After Split",&(top->idef));
#endif
}
