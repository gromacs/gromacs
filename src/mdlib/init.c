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
#include "fatal.h"
#include "symtab.h"
#include "txtdump.h"
#include "splittop.h"
#include "mdatoms.h"
#include "mdrun.h"
#include "statutil.h"
#include "names.h"
#include "calcgrid.h"

#define BUFSIZE	256

#define NOT_FINISHED(l1,l2) \
  printf("not finished yet: lines %d .. %d in %s\n",l1,l2,__FILE__)

void check_nnodes_top(char *fn,t_topology *top)
{
  int i,j;
  bool bOld=FALSE;
  
  /* This routine resets old parallel topologies to single processor
   * topologies.
   */
  
  for(i=0; (i<ebNR); i++) {
    top->blocks[i].multinr[0] = top->blocks[i].nr;
    for(j=1; (j<MAXNODES); j++) {
      bOld = bOld || (top->blocks[i].multinr[j] != 0);
      top->blocks[i].multinr[j] = 0;
    }
  }
  for(i=0; (i<F_NRE); i++) {
    top->idef.il[i].multinr[0] = top->idef.il[i].nr;
    for(j=1; (j<MAXNODES); j++) {
      bOld = bOld || (top->idef.il[i].multinr[j] != 0);
      top->idef.il[i].multinr[j] = 0;
    }
  } 
  if (bOld) {
    fprintf((stdlog != NULL) ? stdlog : stderr,
	    "WARNING:\n"
	    "You are using an old multiprocessor tpr file (%s).\n"
	    "It has been converted back to a single processor topology file\n"
	    "before further processing. In case of problems please make\n"
	    "a new tpr file.\n",
	    fn);
  }
}

static char *int_title(char *title,int nodeid,char buf[], int size)
{
  snprintf(buf,size-1,"%s (%d)",title,nodeid);
  
  return buf;
}

void init_single(FILE *log,t_inputrec *inputrec,
		 char *tpxfile,t_topology *top, 
                 t_state *state,t_mdatoms **mdatoms,
		 t_nsborder *nsb)
{
  int         step;
  real        t;
  
  read_tpx_state(tpxfile,&step,&t,inputrec,state,NULL,top);
  check_nnodes_top(tpxfile,top);

  *mdatoms=atoms2md(log,&top->atoms,inputrec->opts.nFreeze,
		    inputrec->eI,inputrec->delta_t,
		    inputrec->bd_fric,inputrec->opts.tau_t,
		    inputrec->efep!=efepNO,FALSE);
  
  pr_inputrec(log,0,"Input Parameters",inputrec);
  calc_nsb(log,&(top->blocks[ebCGS]),1,0,nsb,0);
  print_nsb(log,"Neighbor Search Blocks",nsb);
}

static void distribute_parallel(t_commrec *cr,int left,int right,char *tpxfile)
{
  int         step;
  real        t;
  t_inputrec  inputrec;
  t_topology  top;
  t_state     state;
  int         npmenodes=0;
  
  init_inputrec(&inputrec);
  read_tpx_state(tpxfile,&step,&t,&inputrec,&state,NULL,&top);
  check_nnodes_top(tpxfile,&top);
  mv_data(cr,left,right,&inputrec,&top,&state);
  done_top(&top);
  done_state(&state);
  done_inputrec(&inputrec);
}

void init_parallel(FILE *log,char *tpxfile,t_commrec *cr,
		   t_inputrec *inputrec,t_topology *top,
		   t_state *state,t_mdatoms **mdatoms,
		   int list)
{
  char buf[256];
  
  if (MASTER(cr)) 
    distribute_parallel(cr,cr->left,cr->right,tpxfile);
    
    /* Read the actual data */
  ld_data(cr,cr->left,cr->right,inputrec,top,state);
  if (cr->nodeid != 0)
    mv_data(cr,cr->left,cr->right,inputrec,top,state);

  /* Make sure the random seeds are different on each node */
  inputrec->ld_seed += cr->nodeid;
  
  /* Printing */
  if (list) {
    if (list&LIST_INPUTREC)
      pr_inputrec(log,0,"parameters of the run",inputrec);
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
  *mdatoms=atoms2md(log,&(top->atoms),inputrec->opts.nFreeze,
		    inputrec->eI,inputrec->delta_t,
		    inputrec->bd_fric,inputrec->opts.tau_t,
		    inputrec->efep!=efepNO,FALSE);
}


