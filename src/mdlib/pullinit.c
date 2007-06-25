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
#include "gmx_fatal.h"
#include "macros.h"
#include "rdgroup.h"
#include "symtab.h"
#include "index.h"
#include "confio.h"
#include "pull.h"
#include "pull_internal.h"
#include "string.h"
#include "network.h"
#include "pbc.h"

static void print_info(FILE *log,t_pull *pull) 
{

  fprintf(log,"\n**************************************************\n"
          "                         PULL INFO                      \n"
          "**************************************************\n");

  switch(pull->ePull) {
  case epullUMBRELLA:
    fprintf(log,"RUN TYPE: Umbrella sampling\n");
    break;
  case epullCONSTRAINT:
    fprintf(log,"RUN TYPE: Constraint\n");
    break;
  default:
    gmx_fatal(FARGS,"No such pull type: %d",pull->ePull);
  }

  if (pull->grp[0].nat == 0) {
    fprintf(log,"REFERENCE TYPE: the origin (0,0,0) --> Absolute Coordinates\n");
  } else if (PULL_CYL(pull)) {
    fprintf(log, "REFERENCE TYPE: center of mass of dynamically made groups\n");
    fprintf(log,"Using dynamic reference groups: r1=%8.3f, r0=%8.3f\n",
	    pull->cyl_r1,pull->cyl_r0);
  } else {
    fprintf(log,"REFERENCE TYPE: center of mass of reference group\n");
  }
}

static void read_whole_index(char *indexfile,char ***grpnames,
                             atom_id ***index, int **ngx,int *totalgrps)
{
  t_block *grps;
  char    **gnames;
  int i,j;

  if(!indexfile)
    gmx_fatal(FARGS,"No index file specified");

  grps = init_index(indexfile,&gnames);
  if(grps->nr==0)
    gmx_fatal(FARGS,"No groups in indexfile\n");

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

static void dd_make_local_pull_group(gmx_domdec_t *dd,
				     t_pullgrp *pg,t_mdatoms *md)
{
  int i,ii;
  gmx_ga2la_t *ga2la=NULL;

  ga2la = dd->ga2la;
  pg->nat_loc = 0;
  for(i=0; i<pg->nat; i++) {
    if (ga2la[pg->ind[i]].cell == 0) {
      ii = ga2la[pg->ind[i]].a;
      if (ii < md->start+md->homenr) {
	/* This is a home atom, add it to the local pull group */
	pg->ind_loc[pg->nat_loc] = ii;
	if (pg->weight)
	  pg->weight_loc[pg->nat_loc] = pg->weight[i];
	pg->nat_loc++;
      }
    }
  }
}

void dd_make_local_pull_groups(gmx_domdec_t *dd,t_pull *pull,t_mdatoms *md)
{
  int g;

  if (pull->eGeom != epullgPOS)
    dd_make_local_pull_group(dd,&pull->grp[0],md);
  for(g=1; g<1+pull->ngrp; g++)
    dd_make_local_pull_group(dd,&pull->grp[g],md);
}

static void init_pull_group_index(FILE *log,t_commrec *cr,
				  t_pullgrp *pg,ivec pulldims,
				  t_mdatoms *md,ivec nFreeze[])
{
  int i,ii,d,start,end,nfreeze,ndim;
  real m,w;
  double wmass,wwmass,buf[3];
  bool bDomDec;
  gmx_ga2la_t *ga2la=NULL;

  bDomDec = DOMAINDECOMP(cr);
  if (bDomDec)
    ga2la = cr->dd->ga2la;

  start = md->start;
  end   = md->homenr + start;

  if (PAR(cr)) {
    snew(pg->ind_loc,pg->nat);
    pg->nat_loc = 0;
    if (pg->weight)
      snew(pg->weight_loc,pg->nat);
  } else {
    pg->nat_loc = pg->nat;
    pg->ind_loc = pg->ind;
    pg->weight_loc = pg->weight;
  }

  nfreeze = 0;
  wmass = 0;
  wwmass = 0;
  for(i=0; i<pg->nat; i++) {
    ii = pg->ind[i];
    if (bDomDec) {
      /* Get the local index of this atom */
      if (ga2la[pg->ind[i]].cell == 0)
	ii = ga2la[pg->ind[i]].a;
      else
	ii = -1;
    }
    if (ii >= start && ii < end) {
      if (PAR(cr) && !bDomDec)
	pg->ind_loc[pg->nat_loc++] = ii;
      if (md->cFREEZE) {
	for(d=0; d<DIM; d++)
	  if (pulldims[d] && nFreeze[md->cFREEZE[ii]][d])
	    nfreeze++;
      }
      m = md->massT[ii];
      if (pg->weight)
	w = pg->weight[i];
      else
	w = 1;
      wmass += w*m;
      wwmass += w*w*m;
    }
  }
  if (PAR(cr)) {
    buf[0]  = wmass;
    buf[1]  = wwmass;
    buf[2]  = nfreeze;
    gmx_sumd(2,buf,cr);
    wmass   = buf[0];
    wwmass  = buf[1];
    nfreeze = (int)(buf[2] + 0.5);
  }

  if (nfreeze == 0) {
    pg->wscale = wmass/wwmass;
    pg->invtm  = 1.0/(pg->wscale*wmass);
  } else {
    ndim = 0;
    for(d=0; d<DIM; d++)
      ndim += pulldims[d]*pg->nat;
    if (nfreeze > 0 && nfreeze < ndim)
      fprintf(log,
	      "\nWARNING: In pull group '%s' some, but not all of the degrees of freedom\n"
	      "         that are subject to pulling are frozen.\n"
	      "         For pulling the whole group will be frozen.\n\n",
	      pg->name);
    pg->wscale = 1.0;
    pg->invtm  = 0.0;
  }
  
  if (bDomDec)
    dd_make_local_pull_group(cr->dd,pg,md);
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
      pgrp->nat = ngx[i];
      snew(pgrp->ind,pgrp->nat);
      for(j=0; j<pgrp->nat; j++)
        pgrp->ind[j] = index[i][j];
      bFound = TRUE;
      if(MASTER(cr))
        fprintf(log,"found group %s: %d elements. First: %d\n",
		pgrp->name,pgrp->nat,pgrp->ind[0]+1);
    }
  }

  if(!bFound)
    gmx_fatal(FARGS,"Can't find group %s in the index file",pgrp->name);

  if (pgrp->nweight > 0 && pgrp->nweight != pgrp->nat)
    gmx_fatal(FARGS,"Number of weights (%d) for pull group '%s' does not match the number of atoms (%d)",pgrp->nweight,pgrp->name,pgrp->nat);
}

void init_pull(FILE *log,int nfile,t_filenm fnm[],t_inputrec *ir,
	       rvec *x,t_mdatoms *md,matrix box,t_commrec *cr) 
{
  int  g,j,m,ii;
  dvec tmp;
  double dist;
  char **grpnames;
  atom_id **index;
  int  *ngx;
  int  totalgrps;    /* total number of groups in the index file */
  char buf[256];
  t_pull *pull;
  t_pullgrp *pgrp;

  pull = &ir->pull;

  /* do we have to do any pulling at all? If not return */
  pull->ePull = epullNO;
  if (opt2bSet("-pi",nfile,fnm) == FALSE)
    return;

  pull->ePBC = ir->ePBC;
  switch (pull->ePBC) {
  case epbcNONE: pull->npbcdim = 0; break;
  case epbcXY:   pull->npbcdim = 2; break;
  default:       pull->npbcdim = 3; break;
  }

  read_pullparams(pull, opt2fn("-pi",nfile,fnm), opt2fn("-po",nfile,fnm));

  if (pull->cyl_r1 > pull->cyl_r0)
    gmx_fatal(FARGS,"pull_r1 > pull_r0");

  if (MASTER(cr)) {
    print_info(log,pull);
  }

  /* read the whole index file */
  read_whole_index(opt2fn("-pn",nfile,fnm),&grpnames,&index,&ngx,&totalgrps);

  if (debug) {
    fprintf(debug,"read_whole_index: %d groups total\n",totalgrps);
    for(g=0; g<totalgrps; g++)
      fprintf(debug,"group %i (%s) %d elements\n",g+1,grpnames[g],ngx[g]);
  }


  /* grab the groups that are specified in the param file */
  for(g=1; g<1+pull->ngrp; g++)
    get_pull_index(log,&pull->grp[g],index,ngx,grpnames,totalgrps,cr);

  if (pull->eGeom == epullgPOS || strcmp(pull->grp[0].name,"") == 0) {
    pull->grp[0].nat = 0;
    pull->grp[0].invtm = 0;
  } else {
    get_pull_index(log,&pull->grp[0],index,ngx,grpnames,totalgrps,cr);
  }

  for(g=0; g<pull->ngrp+1; g++) {
    pgrp = &pull->grp[g];
    if (pgrp->nat > 0) {
      /* Set the indices */
      init_pull_group_index(log,cr,pgrp,pull->dim,md,ir->opts.nFreeze);
      
      /* Set the pbc atom */
      if (pgrp->nat > 1) {
	if (pgrp->pbcatom > 0)
	  pgrp->pbcatom -= 1;
	else
	  pgrp->pbcatom = pgrp->ind[pgrp->nat/2];
      } else {
	pull->grp[g].pbcatom = -1;
      }
      if (MASTER(cr))
	fprintf(log,"Pull group %d pbcatom %d\n",g,pgrp->pbcatom);
    } else {
      pgrp->pbcatom = -1;
    }
  }
  
  /* if we use dynamic reference groups, do some initialising for them */
  if (PULL_CYL(pull)) {
    if (pull->grp[0].nat == 0)
      gmx_fatal(FARGS, "Dynamic reference groups are not support when using absolute reference!\n");
    
    snew(pull->dyna,pull->ngrp+1);
    for(g=1; g<1+pull->ngrp; g++) {
      snew(pull->dyna[g].ind,pull->grp[0].nat);
      snew(pull->dyna[g].weight_loc,pull->grp[0].nat);
    }
  }
  
  /* Determine the initial COMs */
  pull_calc_coms(cr,pull,md,x,NULL,box);
  
  /* Print initial COMs */
  if (MASTER(cr)) {
    for(g=0; g<1+pull->ngrp; g++) {
      if (!(g==0 && PULL_CYL(pull))) {
	pgrp = &pull->grp[g];
	fprintf(log,"Pull inv. mass of group %d: %8.6f\n"
		"Pull initial center of mass: %8.3f %8.3f %8.3f\n",
		g,pgrp->invtm,pgrp->x[XX],pgrp->x[YY],pgrp->x[ZZ]);
	if (PULL_CYL(pull))
	  fprintf(log,"Pull initial cylinder group %d atoms. Weighted mass"
		  "for %d:%8.3f\n"
		  "Pull initial center of mass: %8.3f %8.3f %8.3f\n",
		  pull->dyna[g].nat,g,1.0/pull->dyna[g].invtm,
		  pull->dyna[g].x[XX],
		  pull->dyna[g].x[YY],
		  pull->dyna[g].x[ZZ]);
      }
    }
  }
  
  if (pull->bStart) {
    /* Set the initial distances */
    for(g=1; g<1+pull->ngrp; g++) {
      if (PULL_CYL(pull))
	pull_d_pbc_dx(pull->npbcdim,box,
		      pull->grp[g].x,pull->dyna[g].x,tmp);
      else
	pull_d_pbc_dx(pull->npbcdim,box,
		      pull->grp[g].x,pull->grp[0].x,tmp);
      
      switch (pull->eGeom) {
      case epullgDIST:
      case epullgDIR:
      case epullgCYL:
	dist = 0;
	if (pull->eGeom == epullgDIST) {
	  for(m=0; m<DIM; m++)
	    dist += sqr(pull->dim[m]*tmp[m]);
	  dist = sqrt(dist);
	} else {
	  for(m=0; m<DIM; m++)
	    dist += pull->grp[g].vec[m]*tmp[m];
	}
	/* Add the initial distance to init */
	pull->grp[g].init += dist;
	
	if (MASTER(cr))
	  fprintf(log,"Pull init: %8.5f nm, rate: %8.5f nm/ns\n",
		  pull->grp[g].init,pull->grp[g].rate);
	break;
      case epullgPOS:
	/* Add the initial position of the reference to vec */
	dvec_add(pull->grp[g].vec,tmp,pull->grp[g].vec);
	if (MASTER(cr)) {
	  fprintf(log,"Pull init:");
	  for(m=0; m<DIM; m++) {
	    if (pull->dim[m])
	      fprintf(log," %8.5f",pull->grp[g].vec[m]);
	  }
	  fprintf(log,"\n");
	}
	break;
      }
    }
  }

  /* Only do I/O if we are the MASTER */
  if (MASTER(cr)) {
    if (pull->nstxout > 0) {
      pull->out_x = ffopen(opt2fn("-px",nfile,fnm),"w");
      print_pull_header(pull->out_x,pull);
    }
    if (pull->nstfout > 0) {
      pull->out_f = ffopen(opt2fn("-pf",nfile,fnm),"w");
      print_pull_header(pull->out_f,pull);
    }
  }
}
