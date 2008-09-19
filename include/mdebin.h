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
  int    ie,iconrmsd,ib,isvir,ifvir,ipres,ivir,isurft,ipc,itemp,itc,iu,imu;
  int    ivcos,ivisc;
  int    nE,nEg,nEc,nTC,nU;
  int    *igrp;
} t_mdebin;

extern t_mdebin
*init_mdebin(int fp_ene,
	     const gmx_mtop_t *mtop,
	     const t_inputrec *ir);
/* Initiate MD energy bin and write header to energy file. */

extern void upd_mdebin(t_mdebin *md,FILE *fp_dgdl,
		       bool bSum,
		       real tmass,int step,real time,
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
     
extern void print_ebin_header(FILE *log,int steps,real time,real lamb);

extern void print_ebin(int fp_ene,bool bEne,bool bDR,bool bOR,
		       FILE *log,int step,int nsteps,real time,
		       int mode,bool bCompact,
		       t_mdebin *md,t_fcdata *fcd,
		       gmx_groups_t *groups,t_grpopts *opts);

void 
init_state_energyhistory(t_state *      state);

void
update_state_energyhistory(t_state *    state,
						   t_mdebin *   mdebin);

void
restore_energyhistory_from_state(t_mdebin *  mdebin,
								 t_state *   state);

#endif	/* _mdebin_h */

