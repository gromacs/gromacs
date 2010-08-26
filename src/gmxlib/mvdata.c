/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
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

#define   block_bc(cr,   d) gmx_bcast(     sizeof(d),     &(d),(cr))
#define  nblock_bc(cr,nr,d) { if ((nr) > 0) gmx_bcast((nr)*sizeof((d)[0]), (d),(cr)); }
#define    snew_bc(cr,d,nr) { if (!MASTER(cr)) snew((d),(nr)); }
/* Dirty macro with bAlloc not as an argument */
#define nblock_abc(cr,nr,d) { if (bAlloc) snew((d),(nr)); nblock_bc(cr,(nr),(d)); }

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

static void bc_strings_resinfo(const t_commrec *cr,t_symtab *symtab,
			       int nr,t_resinfo *resinfo)
{
  int  i;
  int  *handle;

  snew(handle,nr);
  if (MASTER(cr)) {
    for(i=0; (i<nr); i++)
      handle[i] = lookup_symtab(symtab,resinfo[i].name);
  }
  nblock_bc(cr,nr,handle);

  if (!MASTER(cr)) {
    for (i=0; (i<nr); i++) 
      resinfo[i].name = get_symtab_handle(symtab,handle[i]);
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
  snew_bc(cr,atoms->resinfo,atoms->nres);
  nblock_bc(cr,atoms->nres,atoms->resinfo);
  bc_strings_resinfo(cr,symtab,atoms->nres,atoms->resinfo);
  /* QMMM requires atomtypes to be known on all nodes as well */
  bc_strings(cr,symtab,atoms->nr,&atoms->atomtype);
  bc_strings(cr,symtab,atoms->nr,&atoms->atomtypeB);
}

static void bc_groups(const t_commrec *cr,t_symtab *symtab,
		      int natoms,gmx_groups_t *groups)
{
  int dummy;
  int g,n;

  bc_grps(cr,groups->grps);
  block_bc(cr,groups->ngrpname);
  bc_strings(cr,symtab,groups->ngrpname,&groups->grpname);
  for(g=0; g<egcNR; g++) {
    if (MASTER(cr)) {
      if (groups->grpnr[g]) {
	n = natoms;
      } else {
	n = 0;
      }
    }
    block_bc(cr,n);
    if (n == 0) {
      groups->grpnr[g] = NULL;
    } else {
      snew_bc(cr,groups->grpnr[g],n);
      nblock_bc(cr,n,groups->grpnr[g]);
    }
  }
  if (debug) fprintf(debug,"after bc_groups\n");
}

void bcast_state_setup(const t_commrec *cr,t_state *state)
{
  block_bc(cr,state->natoms);
  block_bc(cr,state->ngtc);
  block_bc(cr,state->nnhpres);
  block_bc(cr,state->nhchainlength);
  block_bc(cr,state->nrng);
  block_bc(cr,state->nrngi);
  block_bc(cr,state->flags);
}

void bcast_state(const t_commrec *cr,t_state *state,gmx_bool bAlloc)
{
  int i,nnht,nnhtp;

  bcast_state_setup(cr,state);

  nnht = (state->ngtc)*(state->nhchainlength); 
  nnhtp = (state->nnhpres)*(state->nhchainlength); 

  if (MASTER(cr)) {
    bAlloc = FALSE;
  }
  if (bAlloc) {
    state->nalloc = state->natoms;
  }
  for(i=0; i<estNR; i++) {
    if (state->flags & (1<<i)) {
      switch (i) {
      case estLAMBDA:  block_bc(cr,state->lambda); break;
      case estBOX:     block_bc(cr,state->box); break;
      case estBOX_REL: block_bc(cr,state->box_rel); break;
      case estBOXV:    block_bc(cr,state->boxv); break;
      case estPRES_PREV: block_bc(cr,state->pres_prev); break;
      case estSVIR_PREV: block_bc(cr,state->svir_prev); break;
      case estFVIR_PREV: block_bc(cr,state->fvir_prev); break;
      case estNH_XI:   nblock_abc(cr,nnht,state->nosehoover_xi); break;
      case estNH_VXI:  nblock_abc(cr,nnht,state->nosehoover_vxi); break;
      case estNHPRES_XI:   nblock_abc(cr,nnhtp,state->nhpres_xi); break;
      case estNHPRES_VXI:  nblock_abc(cr,nnhtp,state->nhpres_vxi); break;
      case estTC_INT:  nblock_abc(cr,state->ngtc,state->therm_integral); break;
      case estVETA:    block_bc(cr,state->veta); break;
      case estVOL0:    block_bc(cr,state->vol0); break;
      case estX:       nblock_abc(cr,state->natoms,state->x); break;
      case estV:       nblock_abc(cr,state->natoms,state->v); break;
      case estSDX:     nblock_abc(cr,state->natoms,state->sd_X); break;
      case estCGP:     nblock_abc(cr,state->natoms,state->cg_p); break;
	  case estLD_RNG:  if(state->nrngi == 1) nblock_abc(cr,state->nrng,state->ld_rng); break;
	  case estLD_RNGI: if(state->nrngi == 1) nblock_abc(cr,state->nrngi,state->ld_rngi); break;
      case estDISRE_INITF: block_bc(cr,state->hist.disre_initf); break;
      case estDISRE_RM3TAV:
          block_bc(cr,state->hist.ndisrepairs);
          nblock_abc(cr,state->hist.ndisrepairs,state->hist.disre_rm3tav);
          break;
      case estORIRE_INITF: block_bc(cr,state->hist.orire_initf); break;
      case estORIRE_DTAV:
          block_bc(cr,state->hist.norire_Dtav);
          nblock_abc(cr,state->hist.norire_Dtav,state->hist.orire_Dtav);
          break;
      default:
          gmx_fatal(FARGS,
                    "Communication is not implemented for %s in bcast_state",
                    est_names[i]);
      }
    }
  }
}

static void bc_ilists(const t_commrec *cr,t_ilist *ilist)
{
  int ftype;

  /* Here we only communicate the non-zero length ilists */
  if (MASTER(cr)) {
    for(ftype=0; ftype<F_NRE; ftype++) {
      if (ilist[ftype].nr > 0) {
	block_bc(cr,ftype);
	block_bc(cr,ilist[ftype].nr);
	nblock_bc(cr,ilist[ftype].nr,ilist[ftype].iatoms);
      }
    }
    ftype = -1;
    block_bc(cr,ftype);
  } else {
    for(ftype=0; ftype<F_NRE; ftype++) {
      ilist[ftype].nr = 0;
    }
    do {
      block_bc(cr,ftype);
      if (ftype >= 0) {
	block_bc(cr,ilist[ftype].nr);
	snew_bc(cr,ilist[ftype].iatoms,ilist[ftype].nr);
	nblock_bc(cr,ilist[ftype].nr,ilist[ftype].iatoms);
      }
    } while (ftype >= 0);
  }

  if (debug) fprintf(debug,"after bc_ilists\n");
}

static void bc_idef(const t_commrec *cr,t_idef *idef)
{
  block_bc(cr,idef->ntypes);
  block_bc(cr,idef->atnr);
  snew_bc(cr,idef->functype,idef->ntypes);
  snew_bc(cr,idef->iparams,idef->ntypes);
  nblock_bc(cr,idef->ntypes,idef->functype);
  nblock_bc(cr,idef->ntypes,idef->iparams);
  block_bc(cr,idef->fudgeQQ);
  bc_ilists(cr,idef->il);
  block_bc(cr,idef->ilsort);
}

static void bc_cmap(const t_commrec *cr, gmx_cmap_t *cmap_grid)
{
	int i,j,nelem,ngrid;
	
	block_bc(cr,cmap_grid->ngrid);
	block_bc(cr,cmap_grid->grid_spacing);
	
	ngrid = cmap_grid->ngrid;
	nelem = cmap_grid->grid_spacing * cmap_grid->grid_spacing;
	
	if(ngrid>0)
	{
		snew_bc(cr,cmap_grid->cmapdata,ngrid);
		
		for(i=0;i<ngrid;i++)
		{
			snew_bc(cr,cmap_grid->cmapdata[i].cmap,4*nelem);
			nblock_bc(cr,4*nelem,cmap_grid->cmapdata[i].cmap);
		}
	}
}

static void bc_ffparams(const t_commrec *cr,gmx_ffparams_t *ffp)
{
  int i;
  
  block_bc(cr,ffp->ntypes);
  block_bc(cr,ffp->atnr);
  snew_bc(cr,ffp->functype,ffp->ntypes);
  snew_bc(cr,ffp->iparams,ffp->ntypes);
  nblock_bc(cr,ffp->ntypes,ffp->functype);
  nblock_bc(cr,ffp->ntypes,ffp->iparams);
  block_bc(cr,ffp->reppow);
  block_bc(cr,ffp->fudgeQQ);
  bc_cmap(cr,&ffp->cmap_grid);
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
    snew_bc(cr,g->QMbasis,g->ngQM);
    snew_bc(cr,g->QMcharge,g->ngQM);
    snew_bc(cr,g->QMmult,g->ngQM);
    snew_bc(cr,g->bSH,g->ngQM);
    snew_bc(cr,g->CASorbitals,g->ngQM);
    snew_bc(cr,g->CASelectrons,g->ngQM);
    snew_bc(cr,g->SAon,g->ngQM);
    snew_bc(cr,g->SAoff,g->ngQM);
    snew_bc(cr,g->SAsteps,g->ngQM);
    
    if (g->ngQM)
    {
        nblock_bc(cr,g->ngQM,g->QMmethod);
        nblock_bc(cr,g->ngQM,g->QMbasis);
        nblock_bc(cr,g->ngQM,g->QMcharge);
        nblock_bc(cr,g->ngQM,g->QMmult);
        nblock_bc(cr,g->ngQM,g->bSH);
        nblock_bc(cr,g->ngQM,g->CASorbitals);
        nblock_bc(cr,g->ngQM,g->CASelectrons);
        nblock_bc(cr,g->ngQM,g->SAon);
        nblock_bc(cr,g->ngQM,g->SAoff);
        nblock_bc(cr,g->ngQM,g->SAsteps);
        /* end of QMMM stuff */
    }
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
    snew_bc(cr,pgrp->weight,pgrp->nweight);
    nblock_bc(cr,pgrp->nweight,pgrp->weight);
  }
}

static void bc_pull(const t_commrec *cr,t_pull *pull)
{
  int g;

  block_bc(cr,*pull);
  snew_bc(cr,pull->grp,pull->ngrp+1);
  for(g=0; g<pull->ngrp+1; g++)
  {
      bc_pullgrp(cr,&pull->grp[g]);
  }
}

static void bc_inputrec(const t_commrec *cr,t_inputrec *inputrec)
{
  gmx_bool bAlloc=TRUE;
  int i;
  
  block_bc(cr,*inputrec);
  snew_bc(cr,inputrec->flambda,inputrec->n_flambda);
  nblock_bc(cr,inputrec->n_flambda,inputrec->flambda);
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

static void bc_moltype(const t_commrec *cr,t_symtab *symtab,
		       gmx_moltype_t *moltype)
{
  bc_string(cr,symtab,&moltype->name);
  bc_atoms(cr,symtab,&moltype->atoms);
  if (debug) fprintf(debug,"after bc_atoms\n");

  bc_ilists(cr,moltype->ilist);
  bc_block(cr,&moltype->cgs);
  bc_blocka(cr,&moltype->excls);
}

static void bc_molblock(const t_commrec *cr,gmx_molblock_t *molb)
{
  gmx_bool bAlloc=TRUE;
  
  block_bc(cr,molb->type);
  block_bc(cr,molb->nmol);
  block_bc(cr,molb->natoms_mol);
  block_bc(cr,molb->nposres_xA);
  if (molb->nposres_xA > 0) {
    snew_bc(cr,molb->posres_xA,molb->nposres_xA);
    nblock_bc(cr,molb->nposres_xA*DIM,molb->posres_xA[0]);
  }
  block_bc(cr,molb->nposres_xB);
  if (molb->nposres_xB > 0) {
    snew_bc(cr,molb->posres_xB,molb->nposres_xB);
    nblock_bc(cr,molb->nposres_xB*DIM,molb->posres_xB[0]);
  }
  if (debug) fprintf(debug,"after bc_molblock\n");
}

static void bc_atomtypes(const t_commrec *cr, t_atomtypes *atomtypes)
{
  int nr;

  block_bc(cr,atomtypes->nr);

  nr = atomtypes->nr;

  snew_bc(cr,atomtypes->radius,nr);
  snew_bc(cr,atomtypes->vol,nr);
  snew_bc(cr,atomtypes->surftens,nr);
  snew_bc(cr,atomtypes->gb_radius,nr);
  snew_bc(cr,atomtypes->S_hct,nr);

  nblock_bc(cr,nr,atomtypes->radius);
  nblock_bc(cr,nr,atomtypes->vol);
  nblock_bc(cr,nr,atomtypes->surftens);
  nblock_bc(cr,nr,atomtypes->gb_radius);
  nblock_bc(cr,nr,atomtypes->S_hct);
}


void bcast_ir_mtop(const t_commrec *cr,t_inputrec *inputrec,gmx_mtop_t *mtop)
{
  int i; 
  if (debug) fprintf(debug,"in bc_data\n");
  bc_inputrec(cr,inputrec);
  if (debug) fprintf(debug,"after bc_inputrec\n");
  bc_symtab(cr,&mtop->symtab);
  if (debug) fprintf(debug,"after bc_symtab\n");
  bc_string(cr,&mtop->symtab,&mtop->name);
  if (debug) fprintf(debug,"after bc_name\n");

  bc_ffparams(cr,&mtop->ffparams);

  block_bc(cr,mtop->nmoltype);
  snew_bc(cr,mtop->moltype,mtop->nmoltype);
  for(i=0; i<mtop->nmoltype; i++) {
    bc_moltype(cr,&mtop->symtab,&mtop->moltype[i]);
  }

  block_bc(cr,mtop->nmolblock);
  snew_bc(cr,mtop->molblock,mtop->nmolblock);
  for(i=0; i<mtop->nmolblock; i++) {
    bc_molblock(cr,&mtop->molblock[i]);
  }

  block_bc(cr,mtop->natoms);

  bc_atomtypes(cr,&mtop->atomtypes);

  bc_block(cr,&mtop->mols);
  bc_groups(cr,&mtop->symtab,mtop->natoms,&mtop->groups);
}
