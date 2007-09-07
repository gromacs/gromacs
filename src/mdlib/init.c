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

#define BUFSIZE	256

#define NOT_FINISHED(l1,l2) \
  printf("not finished yet: lines %d .. %d in %s\n",l1,l2,__FILE__)

static char *int_title(char *title,int nodeid,char buf[], int size)
{
  snprintf(buf,size-1,"%s (%d)",title,nodeid);
  
  return buf;
}

static void correct_state_entries(t_state *state,t_inputrec *ir)
{
  /* The entries in the state in the tpx file might not correspond
   * with what is needed, so we correct this here.
   */
  state->flags = STATE_HAS_X;
  if (state->x == NULL)
    snew(state->x,state->nalloc);
  if (EI_DYNAMICS(ir->eI)) {
    state->flags = STATE_HAS_V;
    if (state->v == NULL)
      snew(state->v,state->nalloc);
  }
  if (ir->eI == eiSD) {
    state->flags |= STATE_HAS_SDX;
    if (state->sd_X == NULL)
      snew(state->sd_X,state->nalloc);
  }
}

void init_single(FILE *log,t_inputrec *inputrec,
		 char *tpxfile,t_topology *top, 
                 t_state *state)
{
  int         step;
  real        t;
  
  read_tpx_state(tpxfile,&step,&t,inputrec,state,NULL,top);
  correct_state_entries(state,inputrec);
  
  pr_inputrec(log,0,"Input Parameters",inputrec,FALSE);
}

static void distribute_parallel(t_commrec *cr,int left,int right,char *tpxfile,
				t_state *state_p)
{
  int         step;
  real        t;
  t_inputrec  inputrec;
  t_topology  top;
  t_state     *state;
  int         npmenodes=0;
  
  if (state_p) {
    state = state_p;
  } else {
    snew(state,1);
  }

  init_inputrec(&inputrec);
  read_tpx_state(tpxfile,&step,&t,&inputrec,state,NULL,&top);

  mv_data(cr,left,right,&inputrec,&top,state);
  done_top(&top);
  if (state_p == NULL) {
    done_state(state);
    sfree(state);
  }
  done_inputrec(&inputrec);
}

void init_parallel(FILE *log,char *tpxfile,t_commrec *cr,
		   t_inputrec *inputrec,t_topology *top,
		   t_state *state,
		   int list)
{
  char buf[256];
  
  if (MASTER(cr)) 
    distribute_parallel(cr,cr->left,cr->right,tpxfile,
			DOMAINDECOMP(cr) ? state : NULL);
    
    /* Read the actual data */
  ld_data(cr,cr->left,cr->right,inputrec,top,state);
  if (cr->nodeid != 0)
    mv_data(cr,cr->left,cr->right,inputrec,top,state);

  correct_state_entries(state,inputrec);

  if (!EI_TPI(inputrec->eI)) {
    /* Make sure the random seeds are different on each node */
    inputrec->ld_seed += cr->nodeid;
  }
  
  /* Printing */
  if (list) {
    if (list&LIST_INPUTREC)
      pr_inputrec(log,0,"parameters of the run",inputrec,FALSE);
    if (list&LIST_X)
      pr_rvecs(log,0,"box",state->box,DIM);
    if (list&LIST_V)
      pr_rvecs(log,0,"boxv",state->boxv,DIM);
    if (list&LIST_X)
      pr_rvecs(log,0,int_title("x",0,buf,255),state->x,top->atoms.nr);
    if (list&LIST_V)
      pr_rvecs(log,0,int_title("v",0,buf,255),state->v,top->atoms.nr);
    if (list&LIST_TOP)
      pr_top(log,0,int_title("topology",cr->nodeid,buf,255),top,TRUE);
    fflush(log);
  }
}


