/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef _mdebin_h
#define _mdebin_h

#include <stdio.h>

#include "gromacs/fileio/enxio.h"
#include "gromacs/legacyheaders/ebin.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/state.h"

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
    double              delta_t;
    t_ebin             *ebin;
    int                 ie, iconrmsd, ib, ivol, idens, ipv, ienthalpy;
    int                 isvir, ifvir, ipres, ivir, isurft, ipc, itemp, itc, itcb, iu, imu;
    int                 ivcos, ivisc;
    int                 nE, nEg, nEc, nTC, nTCP, nU, nNHC;
    int                *igrp;
    char              **grpnms;
    int                 mde_n, mdeb_n;
    real               *tmp_r;
    rvec               *tmp_v;
    gmx_bool            bConstr;
    gmx_bool            bConstrVir;
    gmx_bool            bTricl;
    gmx_bool            bDynBox;
    gmx_bool            bNHC_trotter;
    gmx_bool            bPrintNHChains;
    gmx_bool            bMTTK;
    gmx_bool            bMu; /* true if dipole is calculated */
    gmx_bool            bDiagPres;
    gmx_bool            bVir;
    gmx_bool            bPress;
    gmx_bool            bSurft;
    int                 f_nre;
    int                 epc;
    real                ref_p;
    int                 etc;
    int                 nCrmsd;
    gmx_bool            bEner[F_NRE];
    gmx_bool            bEInd[egNR];
    char              **print_grpnms;

    FILE               *fp_dhdl; /* the dhdl.xvg output file */
    double             *dE;      /* energy components for dhdl.xvg output */
    t_mde_delta_h_coll *dhc;     /* the delta U components (raw data + histogram) */
    real               *temperatures;
} t_mdebin;


/* delta_h block type enum: the kinds of energies written out. */
enum
{
    dhbtDH   = 0, /* delta H BAR energy difference*/
    dhbtDHDL = 1, /* dH/dlambda derivative */
    dhbtEN,       /* System energy */
    dhbtPV,       /* pV term */
    dhbtEXPANDED, /* expanded ensemble statistics */
    dhbtNR
};



t_mdebin *init_mdebin(ener_file_t       fp_ene,
                      const gmx_mtop_t *mtop,
                      const t_inputrec *ir,
                      FILE             *fp_dhdl);
/* Initiate MD energy bin and write header to energy file. */

FILE *open_dhdl(const char *filename, const t_inputrec *ir,
                const output_env_t oenv);
/* Open the dhdl file for output */

/* update the averaging structures. Called every time
   the energies are evaluated. */
void upd_mdebin(t_mdebin       *md,
                gmx_bool        bDoDHDL,
                gmx_bool        bSum,
                double          time,
                real            tmass,
                gmx_enerdata_t *enerd,
                t_state        *state,
                t_lambda       *fep,
                t_expanded     *expand,
                matrix          lastbox,
                tensor          svir,
                tensor          fvir,
                tensor          vir,
                tensor          pres,
                gmx_ekindata_t *ekind,
                rvec            mu_tot,
                gmx_constr_t    constr);

void upd_mdebin_step(t_mdebin *md);
/* Updates only the step count in md */

void print_ebin_header(FILE *log, gmx_int64_t steps, double time, real lamb);

void print_ebin(ener_file_t fp_ene, gmx_bool bEne, gmx_bool bDR, gmx_bool bOR,
                FILE *log,
                gmx_int64_t step, double time,
                int mode, gmx_bool bCompact,
                t_mdebin *md, t_fcdata *fcd,
                gmx_groups_t *groups, t_grpopts *opts);



/* Between .edr writes, the averages are history dependent,
   and that history needs to be retained in checkpoints.
   These functions set/read the energyhistory_t structure
   that is written to checkpoints in checkpoint.c */

/* Set the energyhistory_t data structure from a mdebin structure */
void update_energyhistory(energyhistory_t * enerhist, t_mdebin * mdebin);

/* Read the energyhistory_t data structure to a mdebin structure*/
void restore_energyhistory_from_state(t_mdebin        * mdebin,
                                      energyhistory_t * enerhist);

#ifdef __cplusplus
}
#endif

#endif  /* _mdebin_h */
