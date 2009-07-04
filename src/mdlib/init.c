/*
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "tpxio.h"
#include "smalloc.h"
#include "vec.h"
#include "main.h"
#include "mvdata.h"
#include "gmx_fatal.h"
#include "symtab.h"
#include "txtdump.h"
#include "mdatoms.h"
#include "mdrun.h"
#include "statutil.h"
#include "names.h"
#include "calcgrid.h"
#include "gmx_random.h"
#include "update.h"
#include "mdebin.h"

#define BUFSIZE	256

#define NOT_FINISHED(l1,l2) \
  printf("not finished yet: lines %d .. %d in %s\n",l1,l2,__FILE__)

static char *int_title(const char *title,int nodeid,char buf[], int size)
{
  sprintf(buf,"%s (%d)",title,nodeid);
  
  return buf;
}

static void set_state_entries(t_state *state,t_inputrec *ir,int nnodes)
{
  /* The entries in the state in the tpx file might not correspond
   * with what is needed, so we correct this here.
   */
  state->flags = 0;
  if (ir->efep != efepNO)
    state->flags |= (1<<estLAMBDA);
  state->flags |= (1<<estX);
  if (state->x == NULL)
    snew(state->x,state->nalloc);
  if (EI_DYNAMICS(ir->eI)) {
    state->flags |= (1<<estV);
    if (state->v == NULL)
      snew(state->v,state->nalloc);
  }
  if (ir->eI == eiSD2) {
    state->flags |= (1<<estSDX);
    if (state->sd_X == NULL) {
      /* sd_X is not stored in the tpx file, so we need to allocate it */
      snew(state->sd_X,state->nalloc);
    }
  }
  if (ir->eI == eiCG) {
    state->flags |= (1<<estCGP);
  }
  if (EI_SD(ir->eI) || ir->eI == eiBD || ir->etc == etcVRESCALE) {
    state->nrng  = gmx_rng_n();
    state->nrngi = 1;
    if (EI_SD(ir->eI) || ir->eI == eiBD) {
      /* This will be correct later with DD */
      state->nrng  *= nnodes;
      state->nrngi *= nnodes;
    }
    state->flags |= ((1<<estLD_RNG) | (1<<estLD_RNGI));
    snew(state->ld_rng, state->nrng);
    snew(state->ld_rngi,state->nrngi);
  } else {
    state->nrng = 0;
  }
  if (ir->ePBC != epbcNONE) {
    state->flags |= (1<<estBOX);
    if (PRESERVE_SHAPE(*ir)) {
      state->flags |= (1<<estBOX_REL);
    }
    if (ir->epc == epcPARRINELLORAHMAN) {
      state->flags |= (1<<estBOXV);
    }
    if (ir->epc != epcNO) {
      state->flags |= (1<<estPRES_PREV);
    }
    if (ir->etc == etcNOSEHOOVER) {
      state->flags |= (1<<estNH_XI);
    }
  }
  if (ir->etc == etcNOSEHOOVER || ir->etc == etcVRESCALE) {
    state->flags |= (1<<estTC_INT);
  }

  init_ekinstate(&state->ekinstate,ir);

  init_energyhistory(&state->enerhist);
}

void init_single(FILE *fplog,t_inputrec *inputrec,
		 const char *tpxfile,gmx_mtop_t *mtop, 
                 t_state *state)
{
  read_tpx_state(tpxfile,inputrec,state,NULL,mtop);
  set_state_entries(state,inputrec,1);

  if (fplog)
    pr_inputrec(fplog,0,"Input Parameters",inputrec,FALSE);
}

void init_parallel(FILE *log,const char *tpxfile,t_commrec *cr,
		   t_inputrec *inputrec,gmx_mtop_t *mtop,
		   t_state *state,
		   int list)
{
  char buf[256];
  
  if (MASTER(cr)) {
    init_inputrec(inputrec);
    read_tpx_state(tpxfile,inputrec,state,NULL,mtop);
    /* When we will be doing domain decomposition with separate PME nodes
     * the rng entries will be too large, we correct for this later.
     */
    set_state_entries(state,inputrec,cr->nnodes);
  }
  bcast_ir_mtop(cr,inputrec,mtop);

  if (inputrec->eI == eiBD || EI_SD(inputrec->eI)) {
    /* Make sure the random seeds are different on each node */
    inputrec->ld_seed += cr->nodeid;
  }
  
  /* Printing */
  if (list!=0 && log!=NULL) 
  {
	  if (list&LIST_INPUTREC)
		  pr_inputrec(log,0,"parameters of the run",inputrec,FALSE);
	  if (list&LIST_X)
		  pr_rvecs(log,0,"box",state->box,DIM);
	  if (list&LIST_X)
		  pr_rvecs(log,0,"box_rel",state->box_rel,DIM);
	  if (list&LIST_V)
		  pr_rvecs(log,0,"boxv",state->boxv,DIM);
	  if (list&LIST_X)
		  pr_rvecs(log,0,int_title("x",0,buf,255),state->x,state->natoms);
	  if (list&LIST_V)
		  pr_rvecs(log,0,int_title("v",0,buf,255),state->v,state->natoms);
	  if (list&LIST_TOP)
		  pr_mtop(log,0,int_title("topology",cr->nodeid,buf,255),mtop,TRUE);
	  fflush(log);
  }
}


