/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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

void set_state_entries(t_state *state,const t_inputrec *ir,int nnodes)
{
  int nnhpres;

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
    if (ir->eI == eiCG)
    {
        state->flags |= (1<<estCGP);
        if (state->cg_p == NULL)
        {
            /* cg_p is not stored in the tpx file, so we need to allocate it */
            snew(state->cg_p,state->nalloc);
        }
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
  state->nnhpres = 0;
  if (ir->ePBC != epbcNONE) {
    state->flags |= (1<<estBOX);
    if (PRESERVE_SHAPE(*ir)) {
      state->flags |= (1<<estBOX_REL);
    }
    if ((ir->epc == epcPARRINELLORAHMAN) || (ir->epc == epcMTTK)) {
      state->flags |= (1<<estBOXV);
    }
    if (ir->epc != epcNO) {
      if (IR_NPT_TROTTER(ir)) {
	state->nnhpres = 1;
	state->flags |= (1<<estNHPRES_XI);
	state->flags |= (1<<estNHPRES_VXI);
	state->flags |= (1<<estSVIR_PREV);
	state->flags |= (1<<estFVIR_PREV);
	state->flags |= (1<<estVETA);
	state->flags |= (1<<estVOL0);
      } else {
	state->flags |= (1<<estPRES_PREV);
      }
    }
  }

  if (ir->etc == etcNOSEHOOVER) {
    state->flags |= (1<<estNH_XI);
    state->flags |= (1<<estNH_VXI);
  }
  
  if (ir->etc == etcVRESCALE) {
    state->flags |= (1<<estTC_INT);
  }
  
  init_gtc_state(state,state->ngtc,state->nnhpres,ir->opts.nhchainlength); /* allocate the space for nose-hoover chains */
  init_ekinstate(&state->ekinstate,ir);

  init_energyhistory(&state->enerhist);
}


void init_parallel(FILE *log, t_commrec *cr, t_inputrec *inputrec,
                   gmx_mtop_t *mtop)
{
  bcast_ir_mtop(cr,inputrec,mtop);


  if (inputrec->eI == eiBD || EI_SD(inputrec->eI)) {
    /* Make sure the random seeds are different on each node */
    inputrec->ld_seed += cr->nodeid;
  }
}


