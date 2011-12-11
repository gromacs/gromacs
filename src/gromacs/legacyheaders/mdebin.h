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

#include "typedefs.h"
#include "sysstuff.h"
#include "ebin.h"
#include "enxio.h"
#include "types/state.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The functions & data structures here determine the content for outputting  
   the .edr file; the file format and actual writing is done with functions
   defined in enxio.h */

/* forward declaration */
typedef struct t_mde_delta_h_coll t_mde_delta_h_coll;

/* This is the collection of energy averages collected during mdrun, and to 
   be written out to the .edr file. */
typedef struct {
  double delta_t;
  t_ebin *ebin;
  int    ie,iconrmsd,ib,ivol,idens,ipv,ienthalpy;
  int    isvir,ifvir,ipres,ivir,isurft,ipc,itemp,itc,itcb,iu,imu;
  int    ivcos,ivisc;
  int    nE,nEg,nEc,nTC,nTCP,nU,nNHC;
  int    *igrp;
  char   **grpnms;
  int    mde_n,mdeb_n;
  real   *tmp_r;
  rvec   *tmp_v;
  gmx_bool	 bConstr;
  gmx_bool   bConstrVir;
  gmx_bool   bTricl;
  gmx_bool   bDynBox;
  gmx_bool   bNHC_trotter;
  gmx_bool   bMTTK;
  gmx_bool   bDiagPres;
  int    f_nre;
  int    epc;
  real   ref_p;
  int	 etc;
  int    nCrmsd;
  gmx_bool   bEner[F_NRE];
  gmx_bool   bEInd[egNR];
  char   **print_grpnms;

  FILE   *fp_dhdl; /* the dhdl.xvg output file */
  gmx_bool dhdl_derivatives; /* whether to write the derivatives to dhdl.xvg */
  t_mde_delta_h_coll *dhc; /* the BAR delta U (raw data + histogram) */
} t_mdebin;

t_mdebin *init_mdebin(ener_file_t fp_ene,
                             const gmx_mtop_t *mtop,
                             const t_inputrec *ir,
                             FILE *fp_dhdl);
/* Initiate MD energy bin and write header to energy file. */

FILE *open_dhdl(const char *filename,const t_inputrec *ir,
		       const output_env_t oenv);
/* Open the dhdl file for output */

/* update the averaging structures. Called every time 
   the energies are evaluated. */
void upd_mdebin(t_mdebin *md, 
                       gmx_bool write_dhdl,
		       gmx_bool bSum,
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

void upd_mdebin_step(t_mdebin *md);
/* Updates only the step count in md */
  
void print_ebin_header(FILE *log,gmx_large_int_t steps,double time,real lamb);

void print_ebin(ener_file_t fp_ene,gmx_bool bEne,gmx_bool bDR,gmx_bool bOR,
		       FILE *log,
		       gmx_large_int_t step,double time,
		       int mode,gmx_bool bCompact,
		       t_mdebin *md,t_fcdata *fcd,
		       gmx_groups_t *groups,t_grpopts *opts);



/* Between .edr writes, the averages are history dependent,
   and that history needs to be retained in checkpoints. 
   These functions set/read the energyhistory_t structure
   that is written to checkpoints in checkpoint.c */

/* Set the energyhistory_t data structure from a mdebin structure */
void update_energyhistory(energyhistory_t * enerhist,t_mdebin * mdebin);

/* Read the energyhistory_t data structure to a mdebin structure*/
void restore_energyhistory_from_state(t_mdebin * mdebin,
                                             energyhistory_t * enerhist);

#ifdef __cplusplus
}
#endif

#endif	/* _mdebin_h */

