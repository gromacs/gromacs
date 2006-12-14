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
#include "block_tx.h"

static char **ld_string(const t_commrec *cr,int src,t_symtab *symtab)
{
  int name;
  
  blockrx(cr,src,name);
  return get_symtab_handle(symtab,name);
}

static int ld_strings(const t_commrec *cr,int src,t_symtab *symtab,char ****nm)
{
  int  i,nr;
  int  *handle;
  char ***NM;

  blockrx(cr,src,nr);
  snew(handle,nr);
  nblockrx(cr,src,nr,handle);

  snew(*nm,nr);
  NM=*nm;
  for (i=0; (i<nr); i++) 
    NM[i]=get_symtab_handle(symtab,handle[i]);
  sfree(handle);

  return nr;
}

static void ld_symtab(const t_commrec *cr,int src,t_symtab *symtab)
{
  int i,nr,len;
  
  blockrx(cr,src,symtab->nr);
  nr=symtab->nr;
  snew(symtab->symbuf,1);
  symtab->symbuf->bufsize=nr;
  snew(symtab->symbuf->buf,nr);
  for (i=0; i<nr; i++)
    {
      blockrx(cr,src,len);
      snew(symtab->symbuf->buf[i],len);
      nblockrx(cr,src,len,symtab->symbuf->buf[i]);
    }
}

static void ld_grps(const t_commrec *cr,int src,t_grps grps[])
{
  int i;
  
  for(i=0; (i<egcNR); i++) {
    blockrx(cr,src,grps[i].nr);
    snew(grps[i].nm_ind,grps[i].nr);
    nblockrx(cr,src,grps[i].nr,grps[i].nm_ind);
  }
  for( ; (i<egcNR); i++) {
    grps[i].nr=0;
    grps[i].nm_ind=NULL;
  }
}

static void ld_atoms(const t_commrec *cr,int src,t_symtab *symtab,t_atoms *atoms)
{
  int atomnr,dummy;

  blockrx(cr,src,atoms->nr);
  snew(atoms->atom,atoms->nr);
  nblockrx(cr,src,atoms->nr,atoms->atom);
  atomnr=ld_strings(cr,src,symtab,&atoms->atomname);
  if (atomnr != atoms->nr)
    gmx_incons("Number of atoms to send around does not match");
  atoms->nres=ld_strings(cr,src,symtab,&atoms->resname);
  atoms->ngrpname=ld_strings(cr,src,symtab,&atoms->grpname);
  /* QMMM requires atomtypes to be known on all nodes as well */
  dummy = ld_strings(cr,src,symtab,&atoms->atomtype);
  ld_grps(cr,src,atoms->grps);
}

void ld_state(const t_commrec *cr,int src,t_state *state)
{
  blockrx(cr,src,state->natoms);
  blockrx(cr,src,state->ngtc);
  blockrx(cr,src,state->box);
  blockrx(cr,src,state->boxv);
  blockrx(cr,src,state->pcoupl_mu);
  nblockrx(cr,src,state->ngtc,state->nosehoover_xi);
  nblockrx(cr,src,state->natoms,state->x);
  nblockrx(cr,src,state->natoms,state->v);
}

static void ld_ilist(const t_commrec *cr,int src,t_ilist *ilist)
{
  blockrx(cr,src,ilist->nr);
  snew(ilist->iatoms,ilist->nr);
  nblockrx(cr,src,ilist->nr,ilist->iatoms);
}

static void ld_idef(const t_commrec *cr,int src,t_idef *idef)
{
  int i;
  
  blockrx(cr,src,idef->ntypes);
  blockrx(cr,src,idef->atnr);
  snew(idef->functype,idef->ntypes);
  snew(idef->iparams,idef->ntypes);
  nblockrx(cr,src,idef->ntypes,idef->functype);
  nblockrx(cr,src,idef->ntypes,idef->iparams);
  for(i=0; (i<F_NRE); i++)
    ld_ilist(cr,src,&idef->il[i]);
}

static void ld_grpopts(const t_commrec *cr,int src,t_grpopts *g)
{
  int i,n;
  
  blockrx(cr,src,g->ngtc);
  blockrx(cr,src,g->ngacc);
  blockrx(cr,src,g->ngfrz);
  blockrx(cr,src,g->ngener);
  snew(g->nrdf,g->ngtc);
  snew(g->tau_t,g->ngtc);
  snew(g->ref_t,g->ngtc);
  snew(g->acc,g->ngacc);
  snew(g->nFreeze,g->ngfrz);
  snew(g->egp_flags,g->ngener*g->ngener);

  nblockrx(cr,src,g->ngtc,g->nrdf);
  nblockrx(cr,src,g->ngtc,g->tau_t);
  nblockrx(cr,src,g->ngtc,g->ref_t);
  nblockrx(cr,src,g->ngacc,g->acc);
  nblockrx(cr,src,g->ngfrz,g->nFreeze);
  nblockrx(cr,src,g->ngener*g->ngener,g->egp_flags);
  snew(g->annealing,g->ngtc);
  snew(g->anneal_npoints,g->ngtc);
  snew(g->anneal_time,g->ngtc);
  snew(g->anneal_temp,g->ngtc);
  nblockrx(cr,src,g->ngtc,g->annealing);
  nblockrx(cr,src,g->ngtc,g->anneal_npoints);
  for(i=0;(i<g->ngtc); i++) {
    n = g->anneal_npoints[i];
    if (n > 0) {
      snew(g->anneal_time[i],n);
      snew(g->anneal_temp[i],n);
      nblockrx(cr,src,n,g->anneal_time[i]);
      nblockrx(cr,src,n,g->anneal_temp[i]);
    }
  }

  /* QMMM stuff, see inputrec */
  blockrx(cr,src,g->ngQM);
  snew(g->QMmethod,g->ngQM);
  nblockrx(cr,src,g->ngQM,g->QMmethod);
  snew(g->QMbasis,g->ngQM);
  nblockrx(cr,src,g->ngQM,g->QMbasis);
  snew(g->QMcharge,g->ngQM);
  nblockrx(cr,src,g->ngQM,g->QMcharge);
  snew(g->QMmult,g->ngQM);
  nblockrx(cr,src,g->ngQM,g->QMmult);
  snew(g->bSH,g->ngQM);
  nblockrx(cr,src,g->ngQM,g->bSH);
  snew(g->CASorbitals,g->ngQM);
  nblockrx(cr,src,g->ngQM,g->CASorbitals);
  snew(g->CASelectrons,g->ngQM);
  nblockrx(cr,src,g->ngQM,g->CASelectrons);
  snew(g->SAon,g->ngQM);
  nblockrx(cr,src,g->ngQM,g->SAon);
  snew(g->SAoff,g->ngQM);
  nblockrx(cr,src,g->ngQM,g->SAoff);
  snew(g->SAsteps,g->ngQM);
  nblockrx(cr,src,g->ngQM,g->SAsteps);
  /* end of QMMM stuff */
}

static void ld_cosines(const t_commrec *cr,int src,t_cosines *cs)
{
  blockrx(cr,src,cs->n);
  snew(cs->a,cs->n);
  snew(cs->phi,cs->n);
  if (cs->n > 0) {
    nblockrx(cr,src,cs->n,cs->a);
    nblockrx(cr,src,cs->n,cs->phi);
  }
}

static void ld_inputrec(const t_commrec *cr,int src,t_inputrec *inputrec)
{
  int i;
  
  blockrx(cr,src,*inputrec);
  ld_grpopts(cr,src,&(inputrec->opts));
  for(i=0; (i<DIM); i++) {
    ld_cosines(cr,src,&(inputrec->ex[i]));
    ld_cosines(cr,src,&(inputrec->et[i]));
  }
}

void ld_data(const t_commrec *cr,int left,int right,t_inputrec *inputrec,
	     t_topology *top,t_state *state)
{
  int i;
  
  ld_inputrec(cr,left,inputrec);
  if (debug) fprintf(stdlog,"after ld_inputrec");
  ld_symtab(cr,left,&top->symtab);
  if (debug) fprintf(stdlog,"after ld_symtab");
  top->name=ld_string(cr,left,&top->symtab);
  if (debug) fprintf(stdlog,"after ld_name");
  ld_atoms(cr,left,&top->symtab,&top->atoms);
  if (debug) fprintf(stdlog,"after ld_atoms");
  ld_idef(cr,left,&top->idef);
  if (debug) fprintf(stdlog,"after ld_idef");
  for (i=0; (i<ebNR); i++) 
    ld_block(cr,left,&top->blocks[i]);
  if (debug) fprintf(stdlog,"after ld_block");
  snew(state->nosehoover_xi,inputrec->opts.ngtc);
  snew(state->x,top->atoms.nr);
  snew(state->v,top->atoms.nr);
  ld_state(cr,left,state);
  if (debug) fprintf(stdlog,"after ld_state");
}

static void mv_grpopts(const t_commrec *cr,int dest,t_grpopts *g)
{
  int i,n;
  
  blocktx(cr,dest,g->ngtc);
  blocktx(cr,dest,g->ngacc);
  blocktx(cr,dest,g->ngfrz);
  blocktx(cr,dest,g->ngener);
  nblocktx(cr,dest,g->ngtc,g->nrdf);
  nblocktx(cr,dest,g->ngtc,g->tau_t);
  nblocktx(cr,dest,g->ngtc,g->ref_t);
  nblocktx(cr,dest,g->ngacc,g->acc);
  nblocktx(cr,dest,g->ngfrz,g->nFreeze);
  nblocktx(cr,dest,g->ngener*g->ngener,g->egp_flags);
  nblocktx(cr,dest,g->ngtc,g->annealing);
  nblocktx(cr,dest,g->ngtc,g->anneal_npoints);
  for(i=0;(i<g->ngtc); i++) {
    n = g->anneal_npoints[i];
    if (n > 0) {
      nblocktx(cr,dest,n,g->anneal_time[i]);
      nblocktx(cr,dest,n,g->anneal_temp[i]);
    }
  }
  /* QMMM stuff, see inputrec.h */
  blocktx(cr,dest,g->ngQM);
  nblocktx(cr,dest,g->ngQM,g->QMmethod);
  nblocktx(cr,dest,g->ngQM,g->QMbasis);
  nblocktx(cr,dest,g->ngQM,g->QMcharge);
  nblocktx(cr,dest,g->ngQM,g->QMmult);
  nblocktx(cr,dest,g->ngQM,g->bSH);
  nblocktx(cr,dest,g->ngQM,g->CASorbitals);
  nblocktx(cr,dest,g->ngQM,g->CASelectrons);
  nblocktx(cr,dest,g->ngQM,g->SAon);
  nblocktx(cr,dest,g->ngQM,g->SAoff);
  nblocktx(cr,dest,g->ngQM,g->SAsteps);
  /* end of QMMM stuff */
}

static void mv_cosines(const t_commrec *cr,int dest,t_cosines *cs)
{
  blocktx(cr,dest,cs->n);
  if (cs->n > 0) {
    nblocktx(cr,dest,cs->n,cs->a);
    nblocktx(cr,dest,cs->n,cs->phi);
  }
}

static void mv_inputrec(const t_commrec *cr,int dest,t_inputrec *inputrec)
{
  int i;
  
  blocktx(cr,dest,*inputrec);
  mv_grpopts(cr,dest,&(inputrec->opts));
  for(i=0; (i<DIM); i++) {
    mv_cosines(cr,dest,&(inputrec->ex[i]));
    mv_cosines(cr,dest,&(inputrec->et[i]));
  }
}

static void mv_string(const t_commrec *cr,int dest,t_symtab *symtab,char **s)
{
  int handle;
  
  handle=lookup_symtab(symtab,s);
  blocktx(cr,dest,handle);
}

static void mv_strings(const t_commrec *cr,int dest,t_symtab *symtab,int nr,
                       char ***nm)
{
  int i;
  int *handle;

  snew(handle,nr);
  for(i=0; (i<nr); i++)
    handle[i]=lookup_symtab(symtab,nm[i]);
  blocktx(cr,dest,nr);
  nblocktx(cr,dest,nr,handle);
  sfree(handle);
}

static void mv_symtab(const t_commrec *cr,int dest,t_symtab *symtab)
{
  int i,nr,len;
  struct symbuf *symbuf;
  
  blocktx(cr,dest,symtab->nr);
  nr=symtab->nr;
  symbuf=symtab->symbuf;
  while (symbuf!=NULL)
    {
      for (i=0; (i<symbuf->bufsize)&&(i<nr); i++)
        {
          len=strlen(symbuf->buf[i])+1;
          blocktx(cr,dest,len);
          nblocktx(cr,dest,len,symbuf->buf[i]);
        }
      nr-=i;
      symbuf=symbuf->next;
    }
  if (nr != 0)
    gmx_incons("Sending strings around the ring");
}

static void mv_grps(const t_commrec *cr,int dest,t_grps grps[])
{
  int i;
  
  for(i=0; (i<egcNR); i++) {
    blocktx(cr,dest,grps[i].nr);
    nblocktx(cr,dest,grps[i].nr,grps[i].nm_ind);
  }
}


static void mv_atoms(const t_commrec *cr,int dest,t_symtab *symtab,t_atoms *atoms)
{
  int nr;
  
  nr=atoms->nr;
  blocktx(cr,dest,nr);
  nblocktx(cr,dest,nr,atoms->atom);
  mv_strings(cr,dest,symtab,atoms->nr,atoms->atomname);
  mv_strings(cr,dest,symtab,atoms->nres,atoms->resname);
  mv_strings(cr,dest,symtab,atoms->ngrpname,atoms->grpname);
  /* QMMM requires atomtypes to be know on all nodes */
  mv_strings(cr,dest,symtab,atoms->nr,atoms->atomtype);

  mv_grps(cr,dest,atoms->grps);
}

void mv_state(const t_commrec *cr,int dest,t_state *state)
{
  blocktx(cr,dest,state->natoms);
  blocktx(cr,dest,state->ngtc);
  blocktx(cr,dest,state->box);
  blocktx(cr,dest,state->boxv);
  blocktx(cr,dest,state->pcoupl_mu);
  nblocktx(cr,dest,state->ngtc,state->nosehoover_xi);
  nblocktx(cr,dest,state->natoms,state->x);
  nblocktx(cr,dest,state->natoms,state->v);
}

static void mv_ilist(const t_commrec *cr,int dest,t_ilist *ilist)
{
  blocktx(cr,dest,ilist->nr);
  nblocktx(cr,dest,ilist->nr,ilist->iatoms);
}

static void mv_idef(const t_commrec *cr,int dest,t_idef *idef)
{
  int i;
  
  blocktx(cr,dest,idef->ntypes);
  blocktx(cr,dest,idef->atnr);
  nblocktx(cr,dest,idef->ntypes,idef->functype);
  nblocktx(cr,dest,idef->ntypes,idef->iparams);
  for(i=0; (i<F_NRE); i++)
    mv_ilist(cr,dest,&idef->il[i]);
}

void mv_data(const t_commrec *cr,int left,int right,t_inputrec *inputrec,
             t_topology *top,t_state *state)
{
  int i;
  
  mv_inputrec(cr,right,inputrec);
  mv_symtab(cr,right,&top->symtab);
  mv_string(cr,right,&top->symtab,top->name);
  mv_atoms(cr,right,&top->symtab,&top->atoms);
  mv_idef(cr,right,&top->idef);
  for (i=0; (i<ebNR); i++) 
    mv_block(cr,right,&top->blocks[i]);
  mv_state(cr,right,state);
}

