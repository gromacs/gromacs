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
static char *SRCID_typedefs_c = "$Id$";

#include "smalloc.h"
#include "assert.h"
#include "symtab.h"
#include <string.h>

void init_block(t_block *block)
{
  int i;

  block->nr    = 0;
  block->nra   = 0;
  snew(block->index,1);
  block->index[0] = 0;
  block->a     = NULL;
  for(i=0; (i<MAXNODES); i++)
    block->multinr[i]=0;
}

void init_atom(t_atoms *at)
{
  int i;
  
  init_block(&(at->excl));
  at->nr       = 0;
  at->nres     = 0;
  at->ngrpname = 0;
  at->atom     = NULL;
  at->resname  = NULL;
  at->atomname = NULL;
  at->grpname  = NULL;
  at->pdbinfo  = NULL;
  for(i=0; (i<egcNR); i++) {
    at->grps[i].nr=0;
    at->grps[i].nm_ind=NULL;
  }
}

void init_top (t_topology *top)
{
  int i;
  
  top->name = NULL;
  open_symtab(&top->symtab);
  init_atom (&(top->atoms));
  for (i=0; (i<ebNR); i++)
    init_block(&(top->blocks[i]));
}

void init_inputrec(t_inputrec *ir)
{
  memset(ir,0,(size_t)sizeof(*ir));
}

void stupid_fill(t_block *grp,int natom,bool bOneIndexGroup)
{
  int i;

  snew(grp->a,natom);
  for(i=0; (i<natom); i++)
    grp->a[i]=i;
  grp->nra=natom;
  
  if (bOneIndexGroup) {
    snew(grp->index,2);
    grp->index[0]=0;
    grp->index[1]=natom;
    grp->nr=1;
  }
  else {
    snew(grp->index,natom+1);
    for(i=0; (i<=natom); i++)
      grp->index[i]=i;
    grp->nr=natom;
  }
  grp->multinr[0] = natom-1;
}

void done_block(t_block *block)
{
  block->nr    = 0;
  block->nra   = 0;
  sfree(block->index);
  sfree(block->a);
}

void done_atom (t_atoms *at)
{
  done_block(&(at->excl));
  at->nr       = 0;
  at->nres     = 0;
  sfree(at->atom);
  sfree(at->resname);
  sfree(at->atomname);
}

void done_top(t_topology *top)
{
  int i;
  
  done_atom (&(top->atoms));
  done_symtab(&(top->symtab));
  for (i=0; (i<ebNR); i++)
    done_block(&(top->blocks[i]));
}

void done_inputrec(t_inputrec *ir)
{
  int m;
  
  for(m=0; (m<DIM); m++) {
    if (ir->ex[m].a)   sfree(ir->ex[m].a);
    if (ir->ex[m].phi) sfree(ir->ex[m].phi);
    if (ir->et[m].a)   sfree(ir->et[m].a);
    if (ir->et[m].phi) sfree(ir->et[m].phi);
  }
  if (ir->opts.nrdf)    sfree(ir->opts.nrdf);
  if (ir->opts.ref_t)   sfree(ir->opts.ref_t);
  if (ir->opts.tau_t)   sfree(ir->opts.tau_t);
  if (ir->opts.acc)     sfree(ir->opts.acc);
  if (ir->opts.nFreeze) sfree(ir->opts.nFreeze);
}

void init_t_atoms(t_atoms *atoms, int natoms, bool bPdbinfo)
{
  atoms->nr=natoms;
  atoms->nres=0;
  atoms->ngrpname=0;
  snew(atoms->atomname,natoms);
  snew(atoms->resname,natoms);
  snew(atoms->atom,natoms);
  snew(atoms->grpname,natoms);
  if (bPdbinfo)
    snew(atoms->pdbinfo,natoms);
  else
    atoms->pdbinfo=NULL;
  init_block(&atoms->excl);
}

void free_t_atoms(t_atoms *atoms)
{
  int i;

  for(i=0; i<atoms->nr; i++) {
    sfree(*atoms->atomname[i]);
    *atoms->atomname[i]=NULL;
  }
  sfree(atoms->atomname);
  sfree(atoms->resname);
  sfree(atoms->atom);
  if (atoms->pdbinfo)
    sfree(atoms->pdbinfo);
  atoms->nr=0; 
  atoms->nres=0;
  done_block(&atoms->excl);
}     

