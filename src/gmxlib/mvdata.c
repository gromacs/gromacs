/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sysstuff.h>
#include <string.h>
#include "typedefs.h"
#include "main.h"
#include "mvdata.h"
#include "network.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "symtab.h"
#include "vec.h"
#include "tgroup.h"

#define  block_bc(cr,   d) gmx_bcast(     sizeof(d),     &(d),(cr))
#define nblock_bc(cr,nr,d) gmx_bcast((nr)*sizeof((d)[0]), (d),(cr))
#define   snew_bc(cr,d,nr) { if (!MASTER(cr)) snew((d),(nr)); }

static void bc_string(const t_commrec *cr,t_symtab *symtab,char ***s)
{
  int handle;
  
  if (MASTER(cr)) {
    handle = lookup_symtab(symtab,*s);
  }
  block_bc(cr,handle);
  if (!MASTER(cr)) {
    *s = get_symtab_handle(symtab,handle);
  }
}

static void bc_strings(const t_commrec *cr,t_symtab *symtab,int nr,char ****nm)
{
  int  i;
  int  *handle;
  char ***NM;

  snew(handle,nr);
  if (MASTER(cr)) {
    NM = *nm;
    for(i=0; (i<nr); i++)
      handle[i] = lookup_symtab(symtab,NM[i]);
  }
  nblock_bc(cr,nr,handle);

  if (!MASTER(cr)) {
    snew_bc(cr,*nm,nr);
    NM = *nm;
    for (i=0; (i<nr); i++) 
      (*nm)[i] = get_symtab_handle(symtab,handle[i]);
  }
  sfree(handle);
}

static void bc_symtab(const t_commrec *cr,t_symtab *symtab)
{
  int i,nr,len;
  t_symbuf *symbuf;

  block_bc(cr,symtab->nr);
  nr = symtab->nr;
  snew_bc(cr,symtab->symbuf,1);
  symbuf = symtab->symbuf;
  symbuf->bufsize = nr;
  snew_bc(cr,symbuf->buf,nr);
  for (i=0; i<nr; i++) {
    if (MASTER(cr))
      len = strlen(symbuf->buf[i]) + 1;
    block_bc(cr,len);
    snew_bc(cr,symbuf->buf[i],len);
    nblock_bc(cr,len,symbuf->buf[i]);
  }
}

static void bc_block(const t_commrec *cr,t_block *block)
{
  block_bc(cr,block->nr);
  snew_bc(cr,block->index,block->nr+1);
  nblock_bc(cr,block->nr+1,block->index);
}

static void bc_blocka(const t_commrec *cr,t_blocka *block)
{
  block_bc(cr,block->nr);
  snew_bc(cr,block->index,block->nr+1);
  nblock_bc(cr,block->nr+1,block->index);
  block_bc(cr,block->nra);
  if (block->nra) {
    snew_bc(cr,block->a,block->nra);
    nblock_bc(cr,block->nra,block->a);
  }
}

static void bc_grps(const t_commrec *cr,t_grps grps[])
{
  int i;
  
  for(i=0; (i<egcNR); i++) {
    block_bc(cr,grps[i].nr);
    snew_bc(cr,grps[i].nm_ind,grps[i].nr);
    nblock_bc(cr,grps[i].nr,grps[i].nm_ind);
  }
}

static void bc_atoms(const t_commrec *cr,t_symtab *symtab,t_atoms *atoms)
{
  int dummy;

  block_bc(cr,atoms->nr);
  snew_bc(cr,atoms->atom,atoms->nr);
  nblock_bc(cr,atoms->nr,atoms->atom);
  bc_strings(cr,symtab,atoms->nr,&atoms->atomname);
  block_bc(cr,atoms->nres);
  bc_strings(cr,symtab,atoms->nres,&atoms->resname);
  block_bc(cr,atoms->ngrpname);
  bc_strings(cr,symtab,atoms->ngrpname,&atoms->grpname);
  /* QMMM requires atomtypes to be known on all nodes as well */
  bc_strings(cr,symtab,atoms->nr,&atoms->atomtype);
  bc_grps(cr,atoms->grps);
}

void bcast_state(const t_commrec *cr,t_state *state,bool bAlloc)
{
  block_bc(cr,state->natoms);
  block_bc(cr,state->ngtc);
  block_bc(cr,state->flags);
  block_bc(cr,state->box);
  block_bc(cr,state->boxv);
  block_bc(cr,state->pcoupl_mu);
  if (bAlloc && !MASTER(cr)) {
    snew(state->nosehoover_xi,state->ngtc);
    state->nalloc = state->natoms;
    snew(state->x,state->nalloc);
    snew(state->v,state->nalloc);
    if (state->flags & STATE_HAS_SDX)
      snew(state->sd_X,state->nalloc);
  }
  nblock_bc(cr,state->ngtc,state->nosehoover_xi);
  nblock_bc(cr,state->natoms,state->x);
  nblock_bc(cr,state->natoms,state->v);
  if (state->flags & STATE_HAS_SDX)
    nblock_bc(cr,state->natoms,state->sd_X);
  if (state->flags & STATE_HAS_CGP)
    nblock_bc(cr,state->natoms,state->cg_p);
}

static void bc_ilist(const t_commrec *cr,t_ilist *ilist)
{
  block_bc(cr,ilist->nr);
  snew_bc(cr,ilist->iatoms,ilist->nr);
  nblock_bc(cr,ilist->nr,ilist->iatoms);
}

static void bc_idef(const t_commrec *cr,t_idef *idef)
{
  int i;
  
  block_bc(cr,idef->ntypes);
  block_bc(cr,idef->atnr);
  snew_bc(cr,idef->functype,idef->ntypes);
  snew_bc(cr,idef->iparams,idef->ntypes);
  nblock_bc(cr,idef->ntypes,idef->functype);
  nblock_bc(cr,idef->ntypes,idef->iparams);
  for(i=0; (i<F_NRE); i++)
    bc_ilist(cr,&idef->il[i]);
}

static void bc_grpopts(const t_commrec *cr,t_grpopts *g)
{
  int i,n;
  
  block_bc(cr,g->ngtc);
  block_bc(cr,g->ngacc);
  block_bc(cr,g->ngfrz);
  block_bc(cr,g->ngener);
  snew_bc(cr,g->nrdf,g->ngtc);
  snew_bc(cr,g->tau_t,g->ngtc);
  snew_bc(cr,g->ref_t,g->ngtc);
  snew_bc(cr,g->acc,g->ngacc);
  snew_bc(cr,g->nFreeze,g->ngfrz);
  snew_bc(cr,g->egp_flags,g->ngener*g->ngener);

  nblock_bc(cr,g->ngtc,g->nrdf);
  nblock_bc(cr,g->ngtc,g->tau_t);
  nblock_bc(cr,g->ngtc,g->ref_t);
  nblock_bc(cr,g->ngacc,g->acc);
  nblock_bc(cr,g->ngfrz,g->nFreeze);
  nblock_bc(cr,g->ngener*g->ngener,g->egp_flags);
  snew_bc(cr,g->annealing,g->ngtc);
  snew_bc(cr,g->anneal_npoints,g->ngtc);
  snew_bc(cr,g->anneal_time,g->ngtc);
  snew_bc(cr,g->anneal_temp,g->ngtc);
  nblock_bc(cr,g->ngtc,g->annealing);
  nblock_bc(cr,g->ngtc,g->anneal_npoints);
  for(i=0;(i<g->ngtc); i++) {
    n = g->anneal_npoints[i];
    if (n > 0) {
      snew_bc(cr,g->anneal_time[i],n);
      snew_bc(cr,g->anneal_temp[i],n);
      nblock_bc(cr,n,g->anneal_time[i]);
      nblock_bc(cr,n,g->anneal_temp[i]);
    }
  }

  /* QMMM stuff, see inputrec */
  block_bc(cr,g->ngQM);
  snew_bc(cr,g->QMmethod,g->ngQM);
  nblock_bc(cr,g->ngQM,g->QMmethod);
  snew_bc(cr,g->QMbasis,g->ngQM);
  nblock_bc(cr,g->ngQM,g->QMbasis);
  snew_bc(cr,g->QMcharge,g->ngQM);
  nblock_bc(cr,g->ngQM,g->QMcharge);
  snew_bc(cr,g->QMmult,g->ngQM);
  nblock_bc(cr,g->ngQM,g->QMmult);
  snew_bc(cr,g->bSH,g->ngQM);
  nblock_bc(cr,g->ngQM,g->bSH);
  snew_bc(cr,g->CASorbitals,g->ngQM);
  nblock_bc(cr,g->ngQM,g->CASorbitals);
  snew_bc(cr,g->CASelectrons,g->ngQM);
  nblock_bc(cr,g->ngQM,g->CASelectrons);
  snew_bc(cr,g->SAon,g->ngQM);
  nblock_bc(cr,g->ngQM,g->SAon);
  snew_bc(cr,g->SAoff,g->ngQM);
  nblock_bc(cr,g->ngQM,g->SAoff);
  snew_bc(cr,g->SAsteps,g->ngQM);
  nblock_bc(cr,g->ngQM,g->SAsteps);
  /* end of QMMM stuff */
}

static void bc_cosines(const t_commrec *cr,t_cosines *cs)
{
  block_bc(cr,cs->n);
  snew_bc(cr,cs->a,cs->n);
  snew_bc(cr,cs->phi,cs->n);
  if (cs->n > 0) {
    nblock_bc(cr,cs->n,cs->a);
    nblock_bc(cr,cs->n,cs->phi);
  }
}

static void bc_pullgrp(const t_commrec *cr,t_pullgrp *pgrp)
{
  block_bc(cr,*pgrp);
  if (pgrp->nat > 0) {
    snew_bc(cr,pgrp->ind,pgrp->nat);
    nblock_bc(cr,pgrp->nat,pgrp->ind);
  }
  if (pgrp->nweight > 0) {
    snew_bc(cr,pgrp->ind,pgrp->nweight);
    nblock_bc(cr,pgrp->nweight,pgrp->weight);
  }
}

static void bc_pull(const t_commrec *cr,t_pull *pull)
{
  int g;

  block_bc(cr,*pull);
  snew_bc(cr,pull->grp,pull->ngrp+1);
  for(g=0; g<pull->ngrp+1; g++)
    bc_pullgrp(cr,&pull->grp[g]);
}

static void bc_inputrec(const t_commrec *cr,t_inputrec *inputrec)
{
  int i;
  
  block_bc(cr,*inputrec);
  bc_grpopts(cr,&(inputrec->opts));
  if (inputrec->ePull != epullNO) {
    snew_bc(cr,inputrec->pull,1);
    bc_pull(cr,inputrec->pull);
  }
  for(i=0; (i<DIM); i++) {
    bc_cosines(cr,&(inputrec->ex[i]));
    bc_cosines(cr,&(inputrec->et[i]));
  }
}

void bcast_ir_top(const t_commrec *cr,t_inputrec *inputrec,t_topology *top)
{
  int i; 
  if (debug) fprintf(debug,"in bc_data\n");
  bc_inputrec(cr,inputrec);
  if (debug) fprintf(debug,"after bc_inputrec\n");
  bc_symtab(cr,&top->symtab);
  if (debug) fprintf(debug,"after bc_symtab\n");
  bc_string(cr,&top->symtab,&top->name);
  if (debug) fprintf(debug,"after bc_name\n");
  bc_atoms(cr,&top->symtab,&top->atoms);
  if (debug) fprintf(debug,"after bc_atoms\n");
  bc_idef(cr,&top->idef);
  if (debug) fprintf(debug,"after bc_idef\n");
  bc_block(cr,&top->cgs);
  bc_block(cr,&top->mols);
  bc_blocka(cr,&top->excls);
  if (debug) fprintf(debug,"after bc_block\n");
}
