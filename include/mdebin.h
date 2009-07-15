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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _mdebin_h
#define _mdebin_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "sysstuff.h"
#include "ebin.h"
#include "enxio.h"
#include "types/state.h"

typedef struct {
  t_ebin *ebin;
  int    ie,iconrmsd,ib,ivol,idens,ipv;
  int    isvir,ifvir,ipres,ivir,isurft,ipc,itemp,itc,iu,imu;
  int    ivcos,ivisc;
  int    nE,nEg,nEc,nTC,nU;
  int    *igrp;
  char   **grpnms;
  real   *tmp_r;
  rvec   *tmp_v;
} t_mdebin;

extern t_mdebin
*init_mdebin(int fp_ene,
	     const gmx_mtop_t *mtop,
	     const t_inputrec *ir);
/* Initiate MD energy bin and write header to energy file. */

FILE *open_dhdl(const char *filename,t_inputrec *ir);
/* Open the dhdl file for output */

extern void upd_mdebin(t_mdebin *md,FILE *fp_dhdl,
		       bool bSum,
		       double time,
		       real tmass,
		       gmx_enerdata_t *enerd,
		       t_state *state,
		       matrix  lastbox,
		       tensor svir,
		       tensor fvir,
		       tensor vir,
		       tensor pres,
		       gmx_ekindata_t *ekind,
		       rvec mu_tot,
		       gmx_constr_t constr);

extern void upd_mdebin_step(t_mdebin *md);
/* Updates only the step count in md */
  
extern void print_ebin_header(FILE *log,gmx_step_t steps,double time,real lamb);

extern void print_ebin(int fp_ene,bool bEne,bool bDR,bool bOR,
		       FILE *log,
		       gmx_step_t step,double time,
		       int mode,bool bCompact,
		       t_mdebin *md,t_fcdata *fcd,
		       gmx_groups_t *groups,t_grpopts *opts);

extern void 
init_energyhistory(energyhistory_t * enerhist);

extern void
update_energyhistory(energyhistory_t * enerhist,t_mdebin * mdebin);

extern void
restore_energyhistory_from_state(t_mdebin * mdebin,energyhistory_t * enerhist);

#endif	/* _mdebin_h */

