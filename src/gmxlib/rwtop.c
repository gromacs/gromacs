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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#include <stdio.h>
#include <assert.h>
#include "typedefs.h"
#include "smalloc.h"
#include "names.h"
#include "symtab.h"
#include "txtdump.h"
#include "recio.h"
#include "binio.h"
#include "rwtop.h"
#include "bondf.h"

#define BUFSIZE		256
#define ALIGN		sizeof(double)

static void wr_string(FILE *fp,t_symtab *symtab,char **s)
{
  int handle;
  
  handle=lookup_symtab(symtab,s);
  blockwrite(fp,handle);
}

static char **rd_string(FILE *fp,t_symtab *symtab)
{
  int name;
  
  blockread(fp,name);
  return get_symtab_handle(symtab,name);
}

static void wr_ilist(FILE *fp,t_ilist *ilist)
{
  blockwrite(fp,ilist->nr);
  nblockwrite(fp,MAXPROC,ilist->multinr);
  nblockwrite(fp,ilist->nr,ilist->iatoms);
}

static void wr_idef(FILE *fp,t_idef *idef)
{
  int i;

  blockwrite(fp,idef->atnr);
  blockwrite(fp,idef->pid);
  blockwrite(fp,idef->ntypes);
  nblockwrite(fp,idef->ntypes,idef->functype);
  for (i=0; i<idef->ntypes; i++)
    iparams_blockio(fp,write,idef->functype[i],idef->iparams[i]);
  for(i=0; (i<F_NRE); i++)
    wr_ilist(fp,&idef->il[i]);
}

static void rd_ilist(FILE *fp,t_ilist *ilist)
{
  blockread(fp,ilist->nr);
  nblockread(fp,MAXPROC,ilist->multinr);
  snew(ilist->iatoms,ilist->nr);
  nblockread(fp,ilist->nr,ilist->iatoms);
}

static void rd_idef(FILE *fp,t_idef *idef)
{
  int i;

  blockread(fp,idef->atnr);
  blockread(fp,idef->pid);
  blockread(fp,idef->ntypes);
  snew(idef->functype,idef->ntypes);
  snew(idef->iparams,idef->ntypes);
  nblockread(fp,idef->ntypes,idef->functype);
  
  for (i=0; i<idef->ntypes; i++) {
#ifdef DEBUG
    fprintf(stderr,"func: %s\n",interaction_function[idef->functype[i]].name);
#endif
    iparams_blockio(fp,read,idef->functype[i],idef->iparams[i]);
  }
  for(i=0; (i<F_NRE); i++)
    rd_ilist(fp,&idef->il[i]);
}

static void rm_idef(t_idef *idef)
{
  int i;
  
  sfree(idef->functype);
  sfree(idef->iparams);
  for(i=0; (i<F_NRE); i++)
    sfree(idef->il[i].iatoms);
}

static void wr_strings(FILE *fp,t_symtab *symtab,int nr,char ***nm)
{
  int i;
  
  blockwrite(fp,nr);
  for (i=0; i<nr; i++) wr_string(fp,symtab,nm[i]);
}

static int rd_strings(FILE *fp,t_symtab *symtab,char ****nm)
{
  int i,nr;

  blockread(fp,nr);
  snew(*nm,nr);
  for (i=0; i<nr; i++) (*nm)[i]=rd_string(fp,symtab);
  return nr;
}

static void wr_block(FILE *fp,t_block *block)
{
  nblockwrite(fp,MAXPROC,block->multinr);
  blockwrite(fp,block->nr);
  nblockwrite(fp,block->nr+1,block->index);
  blockwrite(fp,block->nra);
  nblockwrite(fp,block->nra,block->a);  
}

static void rd_block(FILE *fp,t_block *block)
{
  nblockread(fp,MAXPROC,block->multinr);
  blockread(fp,block->nr);
  snew(block->index,block->nr+1);
  nblockread(fp,block->nr+1,block->index);
  blockread(fp,block->nra);
  snew(block->a,block->nra);
  nblockread(fp,block->nra,block->a);  
}

static void rm_block(t_block *block)
{
  sfree(block->index);
  sfree(block->a);
}

static void wr_atom(FILE *fp,t_symtab *symtab,int nr,t_atom *atom)
{
  int i;
  
  for (i=0; i<nr; i++)
    {
      blockwrite(fp,atom[i].type);
      blockwrite(fp,atom[i].typeB);
      blockwrite(fp,atom[i].ptype);
      blockwrite(fp,atom[i].m);
      blockwrite(fp,atom[i].q);
      blockwrite(fp,atom[i].mB);
      blockwrite(fp,atom[i].qB);
      blockwrite(fp,atom[i].resnr);
      nblockwrite(fp,egcNR,atom[i].grpnr);
    }
}

static void wr_grps(FILE *fp,t_grps grps[],int ngroup)
{
  int i;
  
  blockwrite(fp,ngroup);
  for(i=0; (i<ngroup); i++) {
    blockwrite(fp,grps[i].nr);
    nblockwrite(fp,grps[i].nr,grps[i].nm_ind);
  }
}

static void wr_atoms(FILE *fp,t_symtab *symtab,t_atoms *atoms)
{
  blockwrite(fp,atoms->nr);
  wr_atom(fp,symtab,atoms->nr,atoms->atom);
  wr_strings(fp,symtab,atoms->nr,atoms->atomname);
  wr_strings(fp,symtab,atoms->nres,atoms->resname);
  wr_strings(fp,symtab,atoms->ngrpname,atoms->grpname);
  wr_block(fp,&atoms->excl);
  wr_grps(fp,atoms->grps,egcNR);
}

static int rd_grps(FILE *fp,t_grps grps[])
{
  int i,ngroup;
  
  blockread(fp,ngroup);
  for(i=0; (i<ngroup); i++) {
    blockread(fp,grps[i].nr);
    snew(grps[i].nm_ind,grps[i].nr);
    nblockread(fp,grps[i].nr,grps[i].nm_ind);
  }
  return ngroup;
}

static void rd_atom(FILE *fp,t_symtab *symtab,int nr,t_atom *atom)
{
  int i;
  
  for (i=0; (i<nr); i++) {
    blockread(fp,atom[i].type);
    blockread(fp,atom[i].typeB);
    blockread(fp,atom[i].ptype);
    blockread(fp,atom[i].m);
    blockread(fp,atom[i].q);
    blockread(fp,atom[i].mB);
    blockread(fp,atom[i].qB);
    blockread(fp,atom[i].resnr);
    nblockread(fp,egcNR,atom[i].grpnr);
  }
}

static void rd_atoms(FILE *fp,t_symtab *symtab,t_atoms *atoms)
{
  int atomnr;

  blockread(fp,atoms->nr);
  snew(atoms->atom,atoms->nr);
  rd_atom(fp,symtab,atoms->nr,atoms->atom);
  atomnr=rd_strings(fp,symtab,&atoms->atomname);
  assert(atomnr==atoms->nr);
  atoms->nres=rd_strings(fp,symtab,&atoms->resname);
  atoms->ngrpname=rd_strings(fp,symtab,&atoms->grpname);
  rd_block(fp,&atoms->excl);
  (void) rd_grps(fp,atoms->grps);
}

static void rm_atoms(t_atoms *atoms)
{
  sfree(atoms->atom);
  sfree(atoms->atomname);
  sfree(atoms->resname);
  sfree(atoms->grpname);
  rm_block(&atoms->excl);
}

long wr_top(FILE *fp,t_topology *top)
{
  int i;
  long fpos;

  fpos=ftell(fp);
  wr_symtab(fp,&top->symtab);
  wr_string(fp,&top->symtab,top->name);
  wr_atoms(fp,&top->symtab,&top->atoms);
  for (i=0; i<ebNR; i++) wr_block(fp,&top->blocks[i]);
  wr_idef(fp,&top->idef);
  return (ftell(fp)-fpos);
}

long rd_top(FILE *fp,t_topology *top)
{
  int i;
  long fpos;

  fpos=ftell(fp);
  rd_symtab(fp,&top->symtab);
  top->name=rd_string(fp,&top->symtab);
  rd_atoms(fp,&top->symtab,&top->atoms);
  for (i=0; i<ebNR; i++) rd_block(fp,&top->blocks[i]);
  rd_idef(fp,&top->idef);
  return (ftell(fp)-fpos);
}

void rm_top(t_topology *top)
{
  int i;

  rm_symtab(&top->symtab);
  rm_atoms(&top->atoms);
  for (i=0; i<ebNR; i++) 
    rm_block(&top->blocks[i]);
  rm_idef(&top->idef);
}

