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
 * GROwing Monsters And Cloning Shrimps
 */
/* This file is threadsafe - keep it that way */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "princ.h"
#include <stdlib.h>
#include "readinp.h"
#include "sysstuff.h"
#include "futil.h"
#include "statutil.h"
#include "vec.h"
#include "smalloc.h"
#include "typedefs.h"
#include "names.h"
#include "fatal.h"
#include "macros.h"
#include "rdgroup.h"
#include "symtab.h"
#include "index.h"
#include "confio.h"
#include "pull.h"
#include "pull_internal.h"
#include "string.h"
#include "pbc.h"

static void print_info(FILE *log,t_pull *pull) 
{

  fprintf(log,"\n**************************************************\n"
          "                         PULL INFO                      \n"
          "**************************************************\n");

  switch(pull->runtype) {
  case eAfm:
    fprintf(log,"RUN TYPE: Afm\n");
    break;
  case eConstraint:
    fprintf(log,"RUN TYPE: Constraint\n");
    break;
  case eUmbrella:
    fprintf(log,"RUN TYPE: Umbrella sampling\n");
    break;
  default:
    fprintf(log,"RUN TYPE: WARNING! pullinit does not know this runtype\n");
  }

  switch(pull->reftype) {
  case eCom: 
    fprintf(log,"REFERENCE TYPE: center of mass of reference group\n");
    break;
  case eComT0:
    fprintf(log,
	    "REFERENCE TYPE: center of mass of reference group at t=0\n");
    break;
  case eDyn:
    fprintf(log, "REFERENCE TYPE: center of mass of dynamically made groups\n");
    fprintf(log,"Using dynamic reference groups: r=%8.3f, rc=%8.3f\n",
	    pull->r,pull->rc);
    break;
  case eDynT0:
    fprintf(log,
	    "REFERENCE TYPE: center of mass of dynamically made groups,\n"
	    "based on the positions of its atoms at t=0\n");
    fprintf(log,"Using dynamic reference groups: r=%8.3f, rc=%8.3f\n",
	    pull->r,pull->rc);
    break;
  default:
    fprintf(log,"REFERENCE TYPE: no clue! What did you do wrong?\n");
  }
}

static void read_whole_index(char *indexfile,char ***grpnames,
                             atom_id ***index, int **ngx,int *totalgrps)
{
  t_block *grps;
  char    **gnames;
  int i,j;

  if(!indexfile)
    fatal_error(0,"No index file specified");

  grps = init_index(indexfile,&gnames);
  if(grps->nr==0)
    fatal_error(0,"No groups in indexfile\n");

  *totalgrps = grps->nr;
  snew(*index,grps->nr);
  snew(*grpnames,grps->nr);
  snew((*ngx),grps->nr); 
  /* memory leak, can't free ngx anymore. 4bytes/group */

  for(i=0; (i<*totalgrps); i++) {
    (*grpnames)[i]=strdup(gnames[i]);
    (*ngx)[i]=grps->index[i+1]-grps->index[i];
    snew((*index)[i],(*ngx)[i]);
    for(j=0; (j<(*ngx)[i]); j++)
      (*index)[i][j]=grps->a[grps->index[i]+j];
  }

  for(i=0; (i<grps->nr); i++)
    sfree(gnames[i]);
  sfree(gnames);
  done_block(grps);
  sfree(grps);
}

static void print_whole_index(char **grpnames, atom_id **index, int *ngx, 
                              int ngrps, t_commrec * cr) 
{
  int i,j;
  FILE *tmp;
  
  if(MASTER(cr)) {
    tmp = ffopen("indexcheck","w");
    for(i=0;i<ngrps;i++) {
      fprintf(tmp,"\nGrp %d: %s. %d elements\n",i,grpnames[i],ngx[i]);
      for(j=0;j<ngx[i];j++)
        fprintf(tmp," %d ",index[i][j]);
    }
    fflush(tmp);
  }
}

static void set_mass(FILE *log,t_pullgrp *pg,ivec pulldims,
		     t_mdatoms *md,ivec nFreeze[])
{
  int i,d;
  real m,w;
  double wmass,wwmass;
  bool bFrozen,bPartial;

  bFrozen = FALSE;
  bPartial = FALSE;
  wmass = 0;
  wwmass = 0;
  for(i=0; i<pg->ngx; i++) {
    for(d=0; d<DIM; d++)
      if (pulldims[d]) {
	if (nFreeze[md->cFREEZE[pg->idx[i]]][d])
	  bFrozen = TRUE;
	else if (bFrozen)
	  bPartial = TRUE;
      }
    m = md->massT[pg->idx[i]];
    if (pg->nweight > 0)
      w = pg->weight[i];
    else
      w = 1;
    wmass += w*m;
    wwmass += w*w*m;
  }
  if (bPartial)
    fprintf(log,
	    "\nWARNING: In pull group '%s' some, but not all of the degrees of freedom\n"
	    "         that are subject to pulling are frozen.\n"
	    "         For pulling the whole group will be frozen.\n\n",
	    pg->name);
  if (bFrozen) {
    pg->wscale = 1.0;
    pg->invtm  = 0.0;
  } else {
    pg->wscale = wmass/wwmass;
    pg->invtm  = 1.0/(pg->wscale*wmass);
  }
}

static void get_pull_index(FILE *log,t_pullgrp *pgrp,
			   atom_id **index,int *ngx,char **grpnames,
			   int totalgrps,t_commrec *cr)
{
  int i,j;
  bool bFound = FALSE;

  if (MASTER(cr))
    fprintf(log,"Looking for group %s: ",pgrp->name);

  for(i=0;i<totalgrps;i++) {
    if(strcmp(pgrp->name,grpnames[i]) == 0) {
      /* found the group we're looking for */
      snew(pgrp->idx,ngx[i]);
      for(j=0;j<ngx[i];j++)
        pgrp->idx[j] = index[i][j];
      pgrp->ngx = ngx[i];
      bFound = TRUE;
      if(MASTER(cr))
        fprintf(log,"found group %s: %d elements. First: %d\n",
		pgrp->name,ngx[i],pgrp->idx[0]);
    }
  }

  if(!bFound)
    fatal_error(0,"Can't find group %s in the index file",pgrp->name);

  if (pgrp->nweight > 0 && pgrp->nweight != pgrp->ngx)
    fatal_error(0,"Number of weights (%d) for pull group '%s' does not match the number of atoms (%d)",pgrp->nweight,pgrp->name,pgrp->ngx);
}
void init_pull(FILE *log,int nfile,t_filenm fnm[],t_pull *pull,rvec *x,
               t_mdatoms *md,ivec nFreeze[],matrix box,
	       int start,int homenr,t_commrec *cr) 
{
  int  i,j,m,ii;
  dvec tmp;
  char **grpnames;
  atom_id **index;
  int  *ngx;
  int  totalgrps;    /* total number of groups in the index file */
  char buf[256];

  /* do we have to do any pulling at all? If not return */
  pull->bPull = opt2bSet("-pi",nfile,fnm);
  if (!pull->bPull)
    return;

  read_pullparams(pull, opt2fn("-pi",nfile,fnm), opt2fn("-po",nfile,fnm));

  if (cr->nnodes > 1 && pull->runtype == eConstraint)
    fatal_error(0,"Can not do constraint force calculation in parallel!\n");

  /* Only do I/O if we are the MASTER */
  if (MASTER(cr)) {
    char * filename = opt2fn("-pd", nfile, fnm); 

    pull->out = ffopen(opt2fn("-pd",nfile,fnm),"w");
  }

  pull->AbsoluteRef = FALSE;

  if(pull->reftype == eDyn || pull->reftype == eDynT0)
    pull->bCyl = TRUE;
  else
    pull->bCyl = FALSE;

  if(pull->bCyl && (pull->rc < 0.01 || pull->r < 0.01))
    fatal_error(0,"rc or r is too small or zero.");

  if(MASTER(cr)) {
    print_info(log,pull);
  }

  /* read the whole index file */
  read_whole_index(opt2fn("-pn",nfile,fnm),&grpnames,&index,&ngx,&totalgrps);

  if(MASTER(cr)) {
    if(pull->bVerbose) {
      fprintf(stderr,"read_whole_index: %d groups total\n",totalgrps);
      for(i=0;i<totalgrps;i++) {
        fprintf(stderr,"group %i (%s) %d elements\n",
                i,grpnames[i],ngx[i]);
      }
    }
  }


  /* grab the groups that are specified in the param file */
  for(i=0;i<pull->ngrp;i++)
    get_pull_index(log,&pull->grp[i],index,ngx,grpnames,totalgrps,cr);

  if (!pull->AbsoluteRef) {
    get_pull_index(log,&pull->ref,index,ngx,grpnames,totalgrps,cr);
    
    /* get more memory! Don't we love C? */
    snew(pull->ref.x0,pull->ref.ngx);
    snew(pull->ref.xp,pull->ref.ngx);
    snew(pull->ref.comhist,pull->reflag);
  } else {
    pull->ref.ngx = 0;
  }

  for(i=0;i<pull->ngrp;i++) {
    set_mass(log,&pull->grp[i],pull->dims,md,nFreeze);
    calc_com(&pull->grp[i],x,md,box);
    copy_dvec(pull->grp[i].x_unc,pull->grp[i].x_con);
    copy_dvec(pull->grp[i].x_unc,pull->grp[i].x_ref);
    copy_dvec(pull->grp[i].x_unc,pull->grp[i].spring);
    if(MASTER(cr)) {
      fprintf(log,"Initializing pull groups. Inv. mass of group %d: %8.3f\n"
              "Initial coordinates center of mass: %8.3f %8.3f %8.3f\n",
              i,pull->grp[i].invtm,
	      pull->grp[i].x_ref[XX],
	      pull->grp[i].x_ref[YY],
	      pull->grp[i].x_ref[ZZ]);
    }
  }

  /* initialize the reference group, in all cases */
  if (!pull->AbsoluteRef) {
    set_mass(log,&pull->ref,pull->dims,md,nFreeze);
    calc_com(&pull->ref,x,md,box);
  } else {
    for(i=0; i<DIM; i++) 
      pull->ref.x_unc[i] = 0;
  }
  copy_dvec(pull->ref.x_unc,pull->ref.x_con);
  copy_dvec(pull->ref.x_unc,pull->ref.x_ref);

  for(j=0;j<pull->reflag;j++)
    copy_dvec(pull->ref.x_unc,pull->ref.comhist[j]);

  if(MASTER(cr))
    fprintf(log,"Initializing reference group. Inv. mass: %8.3f\n"
            "Initial coordinates center of mass: %8.3f %8.3f %8.3f\n",
            pull->ref.invtm,
	    pull->ref.x_ref[XX],
	    pull->ref.x_ref[YY],
	    pull->ref.x_ref[ZZ]);

  /* keep the initial coordinates for center of mass at t0 */
  for(j=0;j<pull->ref.ngx;j++) {
    copy_rvec(x[pull->ref.idx[j]],pull->ref.x0[j]);
    copy_rvec(x[pull->ref.idx[j]],pull->ref.xp[j]);
  }

  /* if we use dynamic reference groups, do some initialising for them */
  if(pull->bCyl) {
    if (pull->AbsoluteRef)
      fatal_error(0, "Dynamic reference groups are not support when using absolute reference!\n");

    snew(pull->dyna,pull->ngrp);
    for(i=0;i<pull->ngrp;i++) {
      snew(pull->dyna[i].idx,pull->ref.ngx);    /* more than nessary */
      snew(pull->dyna[i].weight,pull->ref.ngx);
      snew(pull->dyna[i].comhist,pull->reflag);
    }
    make_refgrps(pull,x,md);
    for(i=0; i<pull->ngrp; i++) {
      copy_dvec(pull->dyna[i].x_unc,pull->dyna[i].x_con);
      copy_dvec(pull->dyna[i].x_unc,pull->dyna[i].x_ref);

      /* initialize comhist values for running averages */
      for(j=0;j<pull->reflag;j++)
        copy_dvec(pull->dyna[i].x_unc,pull->dyna[i].comhist[j]);

      if(MASTER(cr)) {
        fprintf(log,"Initializing dynamic groups. %d atoms. Weighted mass"
                "for %d:%8.3f\n"
                "Initial coordinates center of mass: %8.3f %8.3f %8.3f\n",
                pull->dyna[i].ngx,i,1.0/pull->dyna[i].invtm,
		pull->dyna[i].x_ref[XX],
                pull->dyna[i].x_ref[YY],
		pull->dyna[i].x_ref[ZZ]);
      }
    }
  }

  /* set the reference distances and directions, taking into account pbc */
  for(i=0; i<pull->ngrp; i++) {
    if (pull->bDir) {
      copy_dvec(pull->dir,pull->grp[i].dir);
    } else {
      if(pull->bCyl)
	d_pbc_dx(box,pull->grp[i].x_ref,pull->dyna[i].x_ref,tmp);
      else
	d_pbc_dx(box,pull->grp[i].x_ref,pull->ref.x_ref,tmp);
      
      /* select elements of direction vector to use for Afm and Start runs */
      for(m=0; m<DIM; m++)
	tmp[m] *= pull->dims[m];
      dsvmul(1/dnorm(tmp),tmp,pull->grp[i].dir);
    }

    if(MASTER(cr)) {
      if (pull->runtype == eAfm)
	for(m=0; m<pull->ngrp; ++m) {
	  fprintf(log,"\nPull rate: %e nm/ns. Force constant: %e kJ/(mol nm)",
		  pull->grp[m].AfmRate,pull->grp[m].AfmK);
	  
	  fprintf(log,"\nPull direction: %8.3f %8.3f %8.3f\n",
		  pull->grp[i].dir[XX],pull->grp[i].dir[YY],
		  pull->grp[i].dir[ZZ]);
	}
    }
  }

  if(MASTER(cr))
    print_pull_header(pull);
}
